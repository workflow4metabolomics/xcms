#' wl-24-02-2018, Sat: from https://github.com/RuggeroFerrazza/IsotopicLabelling
#' wl-25-02-2018, Sun: commence for debug and test. Notes: Use other R
#'   packages: ecipex, gsubfn, xcms and stringr
#' wl-08-05-2018, Tue: create an remote repository: 
#'   https://github.com/wanchanglin/isolab 
#' wl-23-03-2019, Sat: apply 'styler' to re-format
#' wl-27-03-2019, Wed: call 'ecipex' directly


#' ========================================================================
#' Main function of the package
#'
#' It computes the estimated X abundances of each sample, returning an
#' object of the class \code{labelling}.
#'
#' @param peak_table A data.frame containing the integrated signals for the
#' samples
#' @param compound The chemical formula of the compound of interest
#' @param charge Natural number, denoting the charge state of the target
#' adduct (1,2,3,...). If not provided, it is 1 by default
#' @param labelling Character, either "H" or "C", specifying which is the
#' labelling element
#' @param mass_shift Maximum shift allowed in the mass range
#' @param RT Expected retention time of the compund of interest
#' @param RT_shift Maximum shift allowed in the retention time range
#' @param chrom_width Chromatographic width of the peaks
#' @param initial_abundance Initial estimate for the abundance of the
#' heaviest X isotope (the label)
#'
#' @details
#' \itemize{
#'  \item{peak_table:}{ The first two columns of \code{peak_table} represent
#'  the mass and the retention time of the peaks; the other columns
#'  represent peak intensities for each sample.  The table can be obtained
#'  using the function \code{table_xcms}}
#'  \item{compound:}{ Character vector, where X has to represent the element
#'  with isotopic distribution to be fitted}
#'  \item{initial_abundance:}{ Numeric vector of the same length as the
#'  number of samples, with the initial estimate for the abundance of the
#'  heaviest X isotope (either ^2H or ^13C). If provided, number between 0
#'  and 100. Otherwise, NA}
#' }
#'
#' @return An object of the class \code{labelling}, which is a list
#' containing the results from the fitting procedure:
#'  \item{compound}{ Character vector specifying the chemical formula of the compound of interest,
#'   with X being the element with unknown isotopic distribution (to be fitted)}
#'  \item{Best_estimate}{ Numeric vector representing the best estimated abundance
#'  of the heaviest X isotope (either ^2H or ^13C). Number between 0 and 100}
#'  \item{std_error}{ Numeric vector containing the standard errors of the estimates}
#'  \item{dev_percent}{ The percentage deviations of the fitted theoretical patterns to the provided experimental patterns}
#'  \item{x_scale}{ Vector containing the \emph{m/z} signals of the isotopic patterns}
#'  \item{y_exp}{ Matrix containing normalised experimental patterns,
#'  where for each sample the most intense signal is set to 100}
#'  \item{y_theor}{ Matrix of normalised fitted theoretical pattern (most intense signal set to 100 for each sample)}
#'  \item{warnings}{ Character vector containing possible warnings coming from the fitting procedure}
#'
#'
#' @examples
#'
#' data(xcms_obj)
#' peak_table <- table_xcms(xcms_obj)
#' fitted_abundances <- main_labelling(peak_table,
#'                                     compound="X40H77NO8P",
#'                                     charge=1,
#'                                     labelling="C",
#'                                     mass_shift=0.05,
#'                                     RT=285,
#'                                     RT_shift=20,
#'                                     chrom_width=7,
#'                                     initial_abundance=NA)
#' summary(object=fitted_abundances)
#' plot(x=fitted_abundances, type="patterns", saveplots=FALSE)
#' plot(x=fitted_abundances, type="residuals", saveplots=FALSE)
#' plot(x=fitted_abundances, type="summary", saveplots=FALSE)
#' save_labelling(fitted_abundances)
#' grouped_estimates <- group_labelling(fitted_abundances,
#'                                      groups=factor(c(rep("C12",4),
#'                                                      rep("C13",4))))
#' # Other possible lipid compounds include:
#' # [PC34:1 + H]+. compound="X42H83NO8P", RT=475, chrom_width=10
#' # [TAG47:3 + NH4]+(a minor species). compound="X50H94NO6",
#' #                                    RT=891, chrom_width=7
#' @author Ruggero Ferrazza
#' @keywords manip
#' @export
main_labelling <- function(peak_table, compound, charge=1, labelling,
                           mass_shift, RT, RT_shift, chrom_width,
                           initial_abundance=NA){

  #' Get some useful isotopic information, to be used in the coming
  #' functions
  info <- isotopic_information(compound, charge, labelling)

  #' Extract one pattern for each sample (each column of the peak_table data
  #' frame)
  experimental_patterns <- isotopic_pattern(peak_table, info, mass_shift, 
                                            RT, RT_shift, chrom_width)

  #' For each extracted pattern, find the X isotopic distribution that
  #' better fits the experimental data
  fitted_abundances <- find_abundance(patterns=experimental_patterns, 
                                      info, initial_abundance, charge)

  return(fitted_abundances)
}


#' ========================================================================
#' Get useful isotopic information
#'
#' This function gathers essential isotopic information required by the
#' other functions of the \code{\link{IsotopicLabelling}} package.
#'
#' @param compound Character vector specifying the chemical formula of the
#'   compound of interest, with X being the element with unknown isotopic
#'   distribution (to be fitted)
#' @param charge Natural number, denoting the charge state of the target
#' adduct (1,2,3,...). If not provided, it is 1 by default
#' @param labelling Character, either "H" or "C", specifying the labelling
#' element
#'
#' @return A list with the following elements:
#' \item{compound}{The same as input}
#' \item{target}{Named vector with the exact masses of all the possible
#' isotopologues arising from the labelling isotope.  M+0 is the
#' monoisotopic mass (sum of the masses of the atoms using the lightest
#' isotope for each element, X included); in M+1 one light isotope is
#' replaced by its heaviest counterpart, and so forth}
#' \item{isotopes}{Table containing the natural isotopic abundances of the
#' elements present in compound (numbers between 0 and 1).
#'  The two isotopes of element X are given NA value}
#' \item{nX}{The number of X atoms. In other words, the number of atoms with
#' unknown isotopic distribution}
#' \item{nTOT}{The total number of atoms of the labelling element (either
#' H+X or C+X)}
#'
#' @details The specified compound is not the neutral molecular species of
#'   interest, but the adduct observed by ESI-MS (such as protonated or
#'   sodiated species). In the chemical formula, the element with unknown
#'   abundance should be denoted by X. For example, the proton adduct of TAG
#'   52:2, C55H103O6, should be written X55H103O6 for ^13C labelling
#'   experiments, and C55X102HO6 for ^2H labelling experiments. Note that in
#'   this last case only 102 hydrogen atoms have unknown isotopic
#'   distribution, since the one giving rise to the adduct comes from the
#'   solvent, and is considered to have fixed natural abundance.
#'
#' @export
#'
#' @examples
#' info <- isotopic_information(compound="X40H77NO8P", charge=1, labelling="C")
#' # This is the case for [PC 32:2+H]+ in a ^13C-labelling experiment
#' @author Ruggero Ferrazza
#' @keywords manip
#'
#' ------------------------------------------------------------------- 
#' wl-20-02-2018, Tue: `nistiso` from package `ecipex`
#' A data frame giving the masses, standard isotopic abundances and nucleon
#' numbers of most chemical elements
isotopic_information <- function(compound, charge=1, labelling){

  #' Check that labelling is correct
  if (labelling !="H" & labelling !="C") 
    stop("Check the labelling character: it should be either H or C")


  #' Introduce X in the isotopes data
  isotopes <- nistiso[,1:3]  # wl-20-02-2018, Tue: global variable from `ecipex`

  X_new <- isotopes[which(isotopes$element==labelling),]
  X_new$element <- "X"

  isotopes <- rbind(isotopes, X_new)


  #' Compute the mass difference between heavier and lighter isotope
  mass_diff <- abs(diff(X_new[,"mass"]))/charge


  #' In isotopes data frame keep only the elements of interest
  DF <- strapply(compound,      #' wl-26-03-18: from package 'gsubfn'
                 "([A-Z][a-z]*)(\\d*)",
                 ~ c(..1, if (nchar(..2)) ..2 else 1),
                 simplify = ~ as.data.frame(t(matrix(..1, 2)), stringsAsFactors = FALSE))
  DF[[2]] <- as.numeric(DF[[2]])


  isotopes <- isotopes[isotopes$element %in% DF[,1],]


  #' Compute the exact mass of the species of interest (suppose that X has natural abundance)
  DF[[3]] <- apply(DF, 1, function(x){
    a <- which(x[1]==isotopes$element)
    a <- a[which.max(isotopes$abundance[a])]
    return(isotopes$mass[a])
  })

  exact_mass <- sum(DF[[2]]*DF[[3]])/charge


  #' Set the X abundance to be unknown
  isotopes$abundance[which(isotopes$element == "X")] <- NA
  row.names(isotopes) <- seq(nrow(isotopes))

  #' Extract total number of atoms of the element being labelled
  nTOT <- DF[which(DF[,1]==labelling),2]
  if (length(nTOT)==0) nTOT <- 0

  nX <- DF[which(DF[,1]=="X"),2]
  if (length(nX)==0) nX <- 0

  nTOT <- nTOT + nX


  #' Create the target vector, containing the exact masses of all the
  #' possible isotopic variants arising from X
  #' Lowest mass: 2 mass units below the monoisotopic mass
  #' Highest mass: 2 mass units above the mass corresponding to all X atoms
  #' having been labelled with the heaviest isotope

  target <- round(seq(from=exact_mass-2*mass_diff, by=mass_diff, length=nX+5), 
                  digits=4)
  names(target) <- c("M-2", "M-1", paste("M+", 0:(nX+2), sep=""))

  return(list(compound=compound, isotopes=isotopes, target=target, nX=nX,
              nTOT=nTOT))

}


#' ========================================================================
#' Extract experimental isotopic patterns from a table of MS peaks
#'
#' Function that extracts the experimental isotopic patterns of a specified
#' compound from a data frame containing MS peak intensities or areas.
#'
#' @param peak_table Data frame of experimental MS peak intensities or areas
#'   (one column for each sample), with the first two columns representing
#'   \emph{m/z} and retention times of the peaks
#' @param info Named list containing isotopic information, output of the
#'   \code{\link{isotopic_information}} function
#' @param mass_shift Maximum difference between theoretical and experimental
#'   mass. In other words, the expected mass accuracy
#' @param RT Expected retention time of the compound of interest
#' @param RT_shift Maximum difference between expected and experimental
#'   retention time of the peaks
#' @param chrom_width An estimate of the chromatographic peak width
#'
#' @return A matrix of extracted experimental isotopic patterns (one column
#'   for each sample), with the first two columns representing the exact
#'   \emph{m/z} and the retention times of the peaks
#'
#' @details The table can be obtained from an \code{xcmsSet} object (output
#'   of the \code{xcms} R package) through the \code{\link{table_xcms}}
#'   function.
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' experimental_patterns <- isotopic_pattern(peak_table, info, mass_shift=0.05,
#' RT=285, RT_shift=20, chrom_width=7)
#' }
#'
#' @author Ruggero Ferrazza
#' @seealso \code{\link{table_xcms}} , \code{\link{isotopic_information}}
#' @keywords manip
isotopic_pattern <- function(peak_table, info, mass_shift, RT, RT_shift, 
                             chrom_width){

  tmp_list <- lapply(info$target, function(x){
    ind <- which( (abs(peak_table$mz - x) < mass_shift) & 
                    (peak_table$rt < (RT + RT_shift) ) & 
                    (peak_table$rt > (RT - RT_shift) ) )
    return(data.frame(ind=ind, rt=peak_table[ind,"rt"]))
  })

  #' Extract the retention times of all the peaks
  rt_overall <- sort(unique(unlist(lapply(tmp_list, function(x){x$rt}), use.names=F)))
  rt_grouped <- apply(abs(outer(rt_overall,rt_overall,'-')), 2, function(u) list(rt_overall[u<=chrom_width]))
  rt_grouped <- unique(lapply(rt_grouped, "[[", 1))

  rt_candidates <- sapply(rt_grouped, mean)
  rt_best <- rt_candidates[which.min(abs(rt_candidates - RT))]

  #' Define the matrix where to put the signals
  patterns <- matrix(0, nrow=length(info$target), ncol=(ncol(peak_table)))
  row.names(patterns) <- names(info$target)
  colnames(patterns) <- colnames(peak_table)


  for (i in 1:length(info$target)){

    ind <- which( (abs(peak_table$mz - info$target[i]) < mass_shift) & 
                    (abs(peak_table$rt - rt_best)<= chrom_width) )

    if (length(ind)>=1){
      ind <- ind[which.min(abs(peak_table$rt[ind] - rt_best))]
      patterns[i,] <- as.numeric(peak_table[ind,])
    }
  }

  #' Check that the most intense signals do not come from M-2 or M-1 (which
  #' could arise from the same species as the target, with one more
  #' unsaturation) If so, "delete" the whole experimental pattern

  if (ncol(patterns)>=3){
    max_pos <- apply(patterns[,-c(1,2)], 2, which.max)
    patterns[,c(F,F,max_pos <=2)]  <- 0
  }

  patterns[is.na(patterns)] <- 0

  #' Cut out the first two masses (M-2 and M-1)
  patterns <- patterns[-c(1,2),]

  #' Check that each pattern has at least two signals different from 0, otherwise set all to 0
  ind_pat <- apply(patterns[,-c(1,2)], 2, function(c) sum(c!=0)) > 1
  patterns[,c(F,F,!ind_pat)] <- 0

  #' Add the exact masses to the patterns matrix
  patterns[,"mz"] <- info$target[-c(1,2)]
  patterns[ which(patterns[,"rt"]==0), "rt"] <- NA

  #' Return the obtained patterns
  return(patterns)

}


#' =======================================================================
#' Fit experimental isotopic patterns
#'
#' Function that takes each of the provided experimental MS isotopic
#' patterns, and fits the best theoretical pattern that reproduces it
#' through a weighted non-linear least squares procedure.
#'
#'
#' @param patterns A matrix of experimental isotopic patterns (one column
#'   for each sample), with the first two columns representing \emph{m/z}
#'   and retention time of the corresponding peaks
#' @param info Named list containing isotopic information, output of the
#' \code{\link{isotopic_information}} function
#' @param initial_abundance Either NA, or a numeric vector of length equal
#' to the number of samples, with the initial guesses on the percentage
#' isotopic abundance of the labelling isotope (denoted as X, it can be
#' either ^2H or ^13C). If provided, numbers between 0 and 100
#' @param charge Natural number, denoting the charge state of the target
#' adduct (1,2,3,...). If not provided, it is 1 by default
#'
#' @return An object of class \code{labelling},
#' which is a list containing the results of the fitting procedure:
#' \item{compound}{Character vector specifying the chemical formula of the compound of interest,
#' with X being the element with unknown isotopic distribution (to be fitted)}
#' \item{best_estimate}{Numeric vector of length equal to the number of samples,
#' containing the estimated percentage abundances of the labelling isotope X
#' (either ^2H or ^13C). Numbers between 0 and 100}
#' \item{std_error}{Numeric vector with the standard errors of the estimates,
#' output of the \code{nls} fitting procedure}
#' \item{dev_percent}{Numeric vector with the percentage deviations between best fitted and related experimental patterns}
#' \item{x_scale}{Numeric vector containing the \emph{m/z} values relative to the signals of the experimental patterns}
#' \item{y_exp}{Matrix of normalised experimental isotopic patterns (one column for each sample).
#' The most intense signal of each pattern is set to 100}
#' \item{y_theor}{Matrix of normalised fitted theoretical isotopic patterns (one column for each sample).
#' The most intense signal of each pattern is set to 100}
#' \item{residuals}{Matrix of residuals: each column is the difference between experimental and best fitted theoretical patterns}
#' \item{warnings}{Character vector with possible warnings from the \code{nls} fitting procedure}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fitted_abundances <- find_abundance(patterns, info, initial_abundance=NA, charge=1)
#' }
#'
#' @author Ruggero Ferrazza
#' @seealso \code{\link{isotopic_information}}
find_abundance <- function(patterns, info, initial_abundance=NA, charge=1){

  tmp_results <- list()
  
  #' ========================================================================
  analysis_X <- function(pattern, info, initial_ab=NA, charge=1){

    #' Create a vector of masses
    target <- info$target[-c(1,2)]

    #' Create the list to return in case of errors
    error_list <- list(compound      = info$compound,
                       best_estimate = NA,
                       std_error     = NA,
                       dev_percent   = NA,
                       x_scale       = target,
                       y_exp         = pattern/max(pattern)*100,
                       y_theor       = rep(NA, times = length(pattern)),
                       residuals     = rep(NA, times = length(pattern)),
                       warnings      = "An error occurred")

    if (sum(pattern)==0) return(error_list)

    #' Normalise the experimental pattern (max intensity = 100)
    pattern <- pattern/max(pattern)*100

    #' Find and store the mass of the most intense signal
    mass_max <- target[which.max(pattern)]

    #' First, rough estimate of the X abundance (either 2H or 13C), using
    #' mass_max and the exact mass If the user inserts a first estimate
    #' (initial_abundance), skip this step
    if (is.na(initial_ab)) 
      initial_ab <-  unname(round(((mass_max - target[1])*charge)/info$nX, digits=3))
    if (initial_ab <0) initial_ab <- 0
    if (initial_ab >1) initial_ab <- 1

    #' Fitting procedure to find the best estimate for the X isotopic
    #' abundance. The signals of the pattern are given weights proportional
    #' to the square root of their intensity, so as to give less importance
    #' to noise. Define a function of only one parameter, the abundance,
    #' which will have to be fitted by the nls function

    pattern_fit <- function(abundance) {
      pattern_from_abundance(abundance, info=info, charge=charge)
    }

    fit <- nls(formula   = pattern ~ pattern_fit(abundance),
               start     = list(abundance = initial_ab),
               control   = list(maxiter = 50, tol = 5e-8, warnOnly = T),
               algorithm = "port",
               weights   = sqrt(pattern),
               na.action = na.exclude,
               lower     = 0,
               upper     = 1)

    if (inherits(try(summary(fit), silent=TRUE),"try-error")) return(error_list)

    warnings <- fit$convInfo$stopMessage

    return(list(best_estimate = summary(fit)$coefficients[1]*100,
                std_error     = summary(fit)$coefficients[2]*100,
                dev_percent   = (sqrt(sum((summary(fit)$residuals)^2)/ sum(pattern^2))*100),
                x_scale       = target,
                y_exp         = pattern,
                y_theor       = pattern_fit(summary(fit)$coefficients[1]),
                residuals     = pattern - pattern_fit(summary(fit)$coefficients[1]),
                warnings      = warnings))

  }


  #' ========================================================================
  for (i in 3:ncol(patterns)){
    tmp_results[[i-2]] <- analysis_X(pattern=patterns[,i], info=info,
                                     initial_ab=initial_abundance[i-2]/100,
                                     charge=charge)
  }
  names(tmp_results) <- colnames(patterns[,-c(1,2)])

  #' Create the output list from the results obtained
  best_estimate <- sapply(tmp_results, "[[", "best_estimate")
  std_error     <- sapply(tmp_results, "[[", "std_error")
  dev_percent   <- sapply(tmp_results, "[[", "dev_percent")
  x_scale       <- patterns[,"mz"]
  y_exp         <- sapply(tmp_results, "[[", "y_exp")
  y_theor       <- sapply(tmp_results, "[[", "y_theor")
  residuals     <- sapply(tmp_results, "[[", "residuals")
  warnings      <- sapply(tmp_results, "[[", "warnings")

  results <- list(compound=info$compound, best_estimate=best_estimate,
                  std_error=std_error, dev_percent=dev_percent,
                  x_scale=x_scale, y_exp=y_exp, y_theor=y_theor,
                  residuals=residuals, warnings=warnings)
  class(results) <- "labelling"

  return(results)

}


#' ========================================================================
#' Compute theoretical isotopic patterns
#'
#' Function that computes the theoretical pattern of the specified compound,
#' given the abundance of the labelling isotope X.
#' The function returns a vector of normalised theoretical intensities,
#' where the maximum signal is set to 100;
#' the related masses are in the "target" vector contained in the input list.
#'
#' @param abundance Isotopic abundance of the labelling isotope X (either
#'   ^2H or ^13C); number between 0 and 1
#' @param info Named list containing isotopic information,
#' output of the \code{\link{isotopic_information}} function
#' @param charge Natural number, denoting the charge state of the target
#'   adduct (1,2,3,...). If not provided, it is 1 by default
#'
#' @return A vector representing the normalised isotopic pattern of the
#'   compound of interest, corresponding to the specified isotopic
#'   distribution
#' @export
#'
#'
#' @author Ruggero Ferrazza
#' @references The function makes use of the \code{\link[ecipex]{ecipex}} R package
#' @seealso \code{\link{isotopic_information}}, \code{\link[ecipex]{ecipex}}
#' @keywords manip
#'
#'
pattern_from_abundance <-function(abundance, info, charge=1){

  #' Get the table containing isotopic information
  isotopes <- info$isotopes

  #' Modify the isotopes table with the abundance specified by the user
  isotopes[which(isotopes$element=="X")[1],"abundance"] <- 1 - abundance
  isotopes[which(isotopes$element=="X")[2],"abundance"] <- abundance

  theoret_pattern <- as.matrix(ecipex(info$compound,
                                      isoinfo = isotopes,
                                      limit = 1e-12,
                                      id = FALSE,
                                      sortby = "mass")[[1]])

  #' Correct masses for charge state
  theoret_pattern[,1] <- theoret_pattern[,1]/charge


  #' Group together signals coming from isotopologues with the same nucleon
  #' number, and assign them the proper position
  theoretical_pattern <- unlist(lapply(info$target[-c(1,2)], function(x){
    ind <- which(abs(x - theoret_pattern[,1]) <0.2/charge)
    return(sum(theoret_pattern[ind,2]))
  }))
  #' Normalise the theoretical pattern
  theoretical_pattern <- theoretical_pattern/max(theoretical_pattern)*100

  return(theoretical_pattern)
}


#' =======================================================================
#' Batch process a set of target analytes
#'
#' This function batch processes LC- or GC-MS data, analyzing the isotopic
#' patterns of a set of target analytes specified by the user. As with other
#' functions of the package, it also requires chromatogrpahic information.
#'
#'
#' @param targets A data frame containing information on the target
#'   analytes. Columns 1 to 8 are: "name_compound" (the name of the target
#'   analytes), "compound", "charge", "labelling", "RT", "RT_shift",
#'   "chrom_width", "mass_shift" (the same parameters for
#'   \code{\link{main_labelling}} function). Columns from 9 onwards contain
#'   the initial estimates for the label abundances (one estimate for each
#'   sample). If not known, a single column of NA values can be entered
#' @param groups A factor containing the name of the group of each sample
#'   analysed; The function will calculate summary statistics for the
#'   samples belonging to the same group
#' @param plot_patterns,plot_residuals,plot_results Whether or not to plot
#'   the patterns, the residuals and a summary of the results. If so, pdf
#'   files are created in the working directory
#' @param save_results Whether to save the results of the estimates. If so,
#'   *.csv files are generated
#'
#' @return batch_grouped_estimates A list as long as the number of target
#'   analytes, containing the summary of the fitted results (group
#'   estimates)
#'
#' @export
#'
#' @examples
#' # Get the sample dataset
#' data("xcms_obj")
#'
#' # Convert the MS data set
#' peak_table <- table_xcms(xcms_obj)
#'
#' # Get the example data frame containing target abalytes
#' data("targets")
#'
#'  # Batch process the data
#'  batch_grouped_estimates <- batch_labelling(targets=targets,
#'  groups=factor(c(rep("C12",4), rep("C13",4))),
#'  plot_patterns=FALSE, plot_residuals=FALSE, plot_results=FALSE, save_results=FALSE)
#'
#' @author Ruggero Ferrazza
#'
#' @seealso \link{main_labelling}, \link{group_labelling},
#' \link{save_labelling}, \link{plot.labelling}
#' -------------------------------------------------------------------- 
#' wl-20-02-2018, Tue: add 'peak_table' as argument. 
batch_labelling <- function(peak_table, targets, groups, plot_patterns=T,
                            plot_residuals=F, plot_results=F, 
                            save_results=F){

  #' attach targets data frame
  attach(targets)

  #' Batch process
  batch_grouped_estimates <- list()

  for (i in 1:length(compound)){
    batch_fitted <- main_labelling(peak_table, 
                                   compound=compound[i],
                                   charge=charge[i], 
                                   labelling=labelling[i],
                                   mass_shift=mass_shift[i], 
                                   RT=RT[i],
                                   RT_shift=RT_shift[i],
                                   chrom_width=chrom_width[i],
                                   initial_abundance=as.numeric(targets[i,8:ncol(targets)]))
    #' plot the patterns
    if (plot_patterns) plot(x=batch_fitted, type="patterns", saveplots=T)
    #' plot the residuals
    if (plot_residuals) plot(x=batch_fitted, type="residuals", saveplots=T)
    #' plot the overall results
    if (plot_results) plot(x=batch_fitted, type="summary", saveplots=T)
    #' save the results to a *.csv file
    if (save_results) save_labelling(batch_fitted)

    #' Group the samples and obtain grouped estimates
    batch_grouped_estimates[[i]] <- group_labelling(batch_fitted, groups=groups)
  }
  names(batch_grouped_estimates) <- name_compound
  detach(targets)
  return(batch_grouped_estimates)
}


#' =======================================================================
#' Compute a single estimated isotopic abundance for each sample group
#'
#' This function groups the fitted abundances in order to give a single
#' estimated value for each sample group, with related standard error of the
#' mean that takes into account both the errors relative to each estimate
#' from the fitting procedure, and the variability across samples.
#'
#' @param fitted_abundances Object of class \code{labelling}.  It contains
#' the results of the isotopic pattern analysis @param groups A factor
#' containing the name of the group of each sample analysed; The function
#' will calculate summary statistics for the samples belonging to the same
#' group
#'
#' @return A data frame containing the summary statistics calculated
#' groupwise.  For each row (a group), it details: \item{N}{The number of
#' samples in that group} \item{Mean}{The averaged estimated percentage
#' isotopic abundance of the labelling isotope} \item{SE mean}{The standard
#' error of the mean} \item{t_crit}{The critical value for a 95\% confidence
#' interval of the t distribution with N-1 degrees of freedom} \item{Lower
#' 95\% CI}{The lower 95\% confidence interval value} \item{Upper 95\%
#' CI}{The upper 95\% confidence interval value}
#'
#' @details For each group, the average is simply computed by considering
#' that the obtained individual values are representative of the population.
#'
#' As for the standard deviations, they are obtained using the law of total
#' variance: the overall variance in each group is the sum of two distinct
#' contributions, the first one related to the uncertainties associated in
#' each sample estimate, and the second one arising from the spread of the
#' estimates (biological variability).
#'
#' @export
#'
#' @examples \dontrun{ grouped_estimates <-
#' group_labelling(fitted_abundances, groups=factor(c(rep("CTRL",4),
#' rep("TRTD",4))))
#'}
#'
#' @author Ruggero Ferrazza @keywords manip
#'
#'
group_labelling <- function(fitted_abundances, groups){

  #' Extract the estimated percentage abundances and the std errors of the
  #' fit
  estimates   <- fitted_abundances$best_estimate
  std_err_fit <- fitted_abundances$std_error

  #' Remove NA's from the data
  ind_NA <- which(is.na(estimates) | is.na(std_err_fit))

  if (length(ind_NA) !=0) {
    estimates   <- estimates[-ind_NA]
    std_err_fit <- std_err_fit[-ind_NA]
    groups      <- groups[-ind_NA]
  }

  #' Compute the average for each group
  avg <- tapply(estimates, groups, mean)

  #' Compute the variance due to biological variability
  var_biol <- tapply(estimates, groups, var)

  #' Compute the variance due to single-sample errors
  var_within <- tapply(std_err_fit, groups,
                       function(x){1/length(x)*sum(x^2)})

  #' Compute the total variance
  var_TOT <- var_biol + var_within

  #' Compute the std error of the MEAN
  N        <- tapply(groups,groups,length)
  std_MEAN <- sqrt( var_TOT / N )

  #' Compute the 95% Confidence intervals
  width <- tapply(groups, groups, function(x){ qt(.975, df=length(x)-1) })
  Lower <- avg - width*std_MEAN
  Upper <- avg + width*std_MEAN

  #' Provide the output
  grouped_estimates <- data.frame(N, avg, std_MEAN, width, Lower, Upper);
  names(grouped_estimates) <- c("N", "Mean", "SE mean", "t_crit", 
                                "Lower 95% CI", "Upper 95% CI")

  return(grouped_estimates)
}


#' ========================================================================
#' Plot method for \code{labelling} objects
#'
#' Produces different types of summary plots for a \code{labelling} object.
#'
#' @param x Object of class \code{labelling}
#' @param type The type of output produced. Available options are
#' "patterns", "residuals", "summary"
#' @param saveplots Should the plots be saved as a PDF file?
#' @param ... Other parameters
#'
#' @return One or more plots
#' @details The default (type 'patterns') plot shows, for each sample in the
#'   class \code{labelling} object, the normalized experimental pattern
#'   superimposed to its fitted theoretical pattern. By setting type to
#'   'residuals', the function plots the residuals (the differences between
#'   experimental and best fitted theoretical patterns). Type 'summary'
#'   produces a summary plot showing the estimated percentage abundances
#'   with related standard errors. If \code{saveplot} is TRUE, the plots are
#'   saved to a PDF file in the working directory.
#' @export
#'
#' @examples
#' \dontrun{
#' plot(x=fitted_abundances, type="patterns", saveplots=TRUE)
#' plot(x=fitted_abundances, type="residuals", saveplots=TRUE)
#' plot(x=fitted_abundances, type="summary", saveplots=TRUE)
#' }
#'
#' @author Ruggero Ferrazza
#' @keywords hplot
plot.labelling <- function(x, type="patterns", saveplots=F, ...){
  fitted_abundances <- x
  plot.new()
  old.par <- par(no.readonly = T)

  sample_name <- names(fitted_abundances$best_estimate)

  if (type=="patterns"){
    #' Plot the results
    if (saveplots==T) {
      pdf(paste(fitted_abundances$compound, "_Isotopic_Patterns", ".pdf", sep=""), width=6.5, height=3)
      par(mar=c(3,3,2.5,0.1), mgp=c(2,0.6,0))
    } else par(mar=c(4,4.5,4,0.1), mfrow=c(3,1), ask=T)

    for (k in 1:length(sample_name)){

      plot(fitted_abundances$x_scale, fitted_abundances$y_exp[,k], type="h",
           main=paste(sample_name[k], " ,   Compound: ", fitted_abundances$compound), 
           xlab="Target mass", ylab="Normalised intensity", ylim=c(0,110), cex.main=1, ...)

      text=paste("Fitted X Abundance: (", sprintf("%1.3f", fitted_abundances$best_estimate[k]), "+/-", sprintf("%1.3f", fitted_abundances$std_error[k]), ")\ %")

      mtext(text, cex=0.8)

      points(fitted_abundances$x_scale, fitted_abundances$y_theor[,k], 
             col=2, pch=16, cex=.5)  #' col="red"
      points(fitted_abundances$x_scale[1], -2.5, pch=17, col="blue")

      legend("top", legend=c("Experimental pattern", "Theoretical pattern"), 
             lty=c(1,0), pch=c(0,16), pt.cex=c(0,0.6), col=c(1,2), cex=0.7, 
             horiz=TRUE)
    }
  } else if (type=="residuals") {
    #' Plot the residuals

    if (saveplots==T) {
      pdf(paste(fitted_abundances$compound, "_Residuals", ".pdf", sep=""), width=6.5, height=3)
      par(mar=c(3,3,2.5,0.1), mgp=c(2,0.6,0))
    } else par(mar=c(4,4.5,4,0.1), mfrow=c(3,1), ask=T)


    for (k in 1:length(sample_name)){

      plot(fitted_abundances$x_scale, fitted_abundances$residuals[,k], 
           type="h", main=paste("Residuals for ", sample_name[k]), 
           xlab="Target mass", ylab="Normalised intensity", cex.main=1, ...)
      text=paste("Fitted X Abundance: (", sprintf("%1.3f", fitted_abundances$best_estimate[k]), "+/-", sprintf("%1.3f", fitted_abundances$std_error[k]), ")\ %")

      mtext(text, cex=0.8)

      abline(h=0, col="gray")

    }

  } else if (type=="summary"){

    if (saveplots==T) pdf(paste(fitted_abundances$compound, "_Summary", ".pdf", sep=""), width=10, height=7)
    par(mar=c(6.4,4.1,4.1,2.1))

    plot(fitted_abundances$best_estimate, type="p", pch=16, cex=0.6, xaxt="n", 
         main="Summary of the Estimated Abundances", 
         ylab="Fitted X abundance  [%]", xlab="", cex.main=1.2, ...)
    mtext(paste("Compound: ", fitted_abundances$compound))
    axis(side=1, at=1:length(sample_name), labels=sample_name, las=2, cex.axis=0.7)

    segments(x0  = 1:length(sample_name),
             y0  = fitted_abundances$best_estimate - fitted_abundances$std_error,
             x1  = 1:length(sample_name),
             y1  = fitted_abundances$best_estimate + fitted_abundances$std_error,
             col = "gray")

  }

  if (saveplots==T) dev.off() else par(old.par)
}


#' ========================================================================
#' Export to csv
#'
#' Function that saves the obtained results to a csv file.
#'
#' @param fitted_abundances Object of class \code{labelling}
#' @param path The directory where to save the csv file.
#' If not specified, the results are saved in the working directory
#'
#' @return The "COMPOUND_Estimated_Abundances.csv" file, containing the
#' results of the analysis.
#' For each sample (one for each row) there are four columns:
#' \enumerate{
#' \item The estimated percentage abundance of the labelling isotope (either
#' ^2H or ^13C);
#' \item The related standard error coming from the \code{nls} fitting
#' procedure;
#' \item The percentage deviation between theoretical and experimental
#' isotopic patterns;
#' \item The outcome message from the fitting procedure, to undersand
#' whether there have been any convergence problems.
#' }
#' @export
#' @author Ruggero Ferrazza
#'
#'
#' @keywords IO
#' @seealso \link{main_labelling}
#'
save_labelling <-function(fitted_abundances, path=getwd()){
  table <- cbind(fitted_abundances$best_estimate,
                 fitted_abundances$std_error,
                 fitted_abundances$dev_percent,
                 fitted_abundances$warnings)
  colnames(table) <- c("Best estimate [%]", "Standard Error [%]",
                       "Percentage deviation [%]", 
                       "Fitting outcome messages/Warnings")

  filename <- paste(fitted_abundances$compound, "_Estimated_Abundances.csv", sep="")
  file <- paste(path, filename, sep="/")
  write.csv(x = table, file=file)
}


#' ========================================================================
#' Summary of the labelling fitting
#'
#' Function that produces a summary of the results from an object of class
#' \code{labelling}.
#'
#' @param object Object of class \code{labelling}, output of either
#' \code{\link{main_labelling}} or \code{\link{find_abundance}} functions
#' @param ... Additional parameters
#'
#' @return \item{results}{Matrix containing a summary of the fitted results.
#'   It has two rows, the first containing the estimated percentage isotopic
#'   abundances of the labelling isotope X (^2H or ^13C), and the second one
#'   containing the standard errors from the fitting procedure}
#' @export
#'
#' @author Ruggero Ferrazza
#'
#'
#' @keywords manip
#'
summary.labelling <- function(object, ...){

  fitted_abundances  <- object
  best_est           <- fitted_abundances$best_estimate
  std_error          <- fitted_abundances$std_error
  results            <- rbind(best_est, std_error)
  row.names(results) <- c("Best Estimate [%]", "Standard Error [%]")

  return(results)
}


#' ======================================================================
#' Process \code{xcmsSet}
#'
#' Function that properly converts an \code{xcmsSet} object,
#' from package \code{xcms}, into a table of peaks.
#'
#' @param xcms_obj An xcmsSet object
#'
#' @return \item{peak_table}{Data frame extracted from the \code{xcmsSet}
#' object.  The first two columns represent mass and retention time of the
#' related peaks}
#'
#' @note The output data frame, required by other functions of the
#' \code{\link{IsotopicLabelling}} R package, can be obtained in a number of
#' other independent ways, such as through proprietary software of the
#' vendor of the MS instrument.
#'
#' @author Ruggero Ferrazza
#'
#' @examples
#' data(xcms_obj)
#' peak_table <- table_xcms(xcms_obj)
#'
#' @keywords manip
#' @export
#' -----------------------------------------------------------------------
#' wl-03-04-2018, Tue: should provide dot arguments for 'peakTable'.
#' Actually this dot argument is for 'groupval'
#' -----------------------------------------------------------------------
table_xcms <- function(xcms_obj){

  #' Check that the file in input is an xcmsSet object
  if (class(xcms_obj) != "xcmsSet") 
    stop("ERROR: The provided object is not an xcmsSet object")

  peak_table <- peakTable(xcms_obj)
  n_files    <- length(sampnames(xcms_obj))
  peak_table <- peak_table[,c(which(colnames(peak_table)=="mz" | colnames(peak_table)=="rt"), (ncol(peak_table)-n_files+1):ncol(peak_table))]

  return(peak_table)
}


#' ======================================================================
#' @title Example data frame for batch processing
#'
#' @description A data frame containing information on target analytes to
#' batch process. The first 8 columns contain the same parameters to input
#' when using \code{\link{main_labelling}} function. See the details in
#' \code{\link{batch_labelling}}
#' @author Ruggero Ferrazza
#' @examples
#' data(targets)
#' "targets"


#' ======================================================================
#' @title Example data set
#'
#' @description An "xcmsSet" object containing the results of the LC-MS
#' analysis of 8 lipid extract obtained in a ^13C isotopic labelling
#' experiment.
#' The first 4 samples were obtained from cell cultures grown under normal
#' conditions (natural ^13C abundance), whereas in the remaining 4 the cells
#' were grown using uniformly-labelled 99\% ^13C glucose.  @format an
#' \code{xcmsSet} object from the package \code{xcms}
#' @source Department of Biochemistry, University of Cambridge - UK
#' @author Dr. Julian L Griffin and Dr. Nyasha Munjoma
#' @examples
#' data("xcms_obj")
#' "xcms_obj"


#' ========================================================================
#' IsotopicLabelling-package
#'
#' The \code{IsotopicLabelling} package allows to analyse the isotopic
#' patterns in MS data obtained in isotopic labelling experiments. From the
#' experimental patterns, the package estimates the isotopic abundance of
#' the stable isotope employed in the labelling experiment (either ^2H or
#' ^13C) inside a specified compound.
#'
#'
#' @section Details: Given a data frame of LC-MS or GC-MS peak intensities
#'   or areas (one column for each sample to analyse), the
#'   \code{\link{IsotopicLabelling}} package first extracts the isotopic
#'   patterns of the specified compound, and then performs an isotopic
#'   pattern analysis to estimate the isotopic abundance of the labelling
#'   isotope. This is performed through a weighted non-linear least squares
#'   fitting procedure, where the resulting estimate is the value for which
#'   the theoretical pattern best reproduces the experimental one. During
#'   the fitting, the experimental signals are given weights proportional to
#'   the square root of their intensity, to correct for the non uniform
#'   variance at different intensity levels. The theoretical patterns are
#'   computed using the \code{\link[ecipex]{ecipex}} R package.
#'
#' @section Block diagram:
#' The isotopic pattern analysis can be divided into the following steps:
#' \enumerate{
#' \item Starting from a class \code{xcmsSet} object (from the \code{xcms} R
#' package), generate a data frame of peak signal intensities or areas, with
#' each column corresponding to a sample.  This step can be avoided if the
#' data frame is already available (obtained by other means);
#' \item Extract from the data frame the experimental isotopic patterns of
#' the specified compound (one pattern for each sample).  In the chemical
#' formula of the compound, the element whose abundance is unknown is called
#' "X";
#' \item Normalise the patterns and estimate the abundance of the label
#' through a weighted non-linear least squares fitting procedure.
#' \item Summarize the results.
#' }
#'
#' @docType package
#' @name IsotopicLabelling
#'
#' @import xcms
#' @import ecipex
#' @import stringr
#' @import gsubfn
#'
#' @author Ruggero Ferrazza, Pietro Franceschi
NULL


#' ========================================================================
#' WL: Slightly modification of function 'plot.labelling' for Galaxy only.
#' 
#' wl-22-05-2018, Tue: For details, see 'plot.labelling'
#' -----------------------------------------------------------------------
plot.func <- function(x, type="patterns", ...){

  sample_name <- names(x$best_estimate)

  if (type=="patterns"){
    old.par <- par(mar=c(4,4.5,4,0.1), mfrow=c(3,1))
    for (k in 1:length(sample_name)){
      plot(x$x_scale, x$y_exp[,k], type="h",
           main=paste(sample_name[k], " ,   Compound: ", x$compound), 
           xlab="Target mass", ylab="Normalised intensity", ylim=c(0,110), cex.main=1, ...)
      text=paste("Fitted X Abundance: (", sprintf("%1.3f", x$best_estimate[k]), "+/-", sprintf("%1.3f", x$std_error[k]), ")\ %")
      mtext(text, cex=0.8)
      points(x$x_scale, x$y_theor[,k], 
             col=2, pch=16, cex=.5)  #' col="red"
      points(x$x_scale[1], -2.5, pch=17, col="blue")
      legend("top", legend=c("Experimental pattern", "Theoretical pattern"), 
             lty=c(1,0), pch=c(0,16), pt.cex=c(0,0.6), col=c(1,2), cex=0.7, 
             horiz=TRUE)
    }
  } else if (type=="residuals") {
    old.par <- par(mar=c(4,4.5,4,0.1), mfrow=c(3,1))
    for (k in 1:length(sample_name)){
      plot(x$x_scale, x$residuals[,k], 
           type="h", main=paste("Residuals for ", sample_name[k]), 
           xlab="Target mass", ylab="Normalised intensity", cex.main=1, ...)
      text=paste("Fitted X Abundance: (", sprintf("%1.3f", x$best_estimate[k]), "+/-", sprintf("%1.3f", x$std_error[k]), ")\ %")
      mtext(text, cex=0.8)
      abline(h=0, col="gray")
    }
  } else if (type=="summary"){
    old.par <- par(mar=c(6.4,4.1,4.1,2.1))
    plot(x$best_estimate, type="p", pch=16, cex=0.6, xaxt="n", 
         main="Summary of the Estimated Abundances", 
         ylab="Fitted X abundance  [%]", xlab="", cex.main=1.2, ...)
    mtext(paste("Compound: ", x$compound))
    axis(side=1, at=1:length(sample_name), labels=sample_name, las=2, cex.axis=0.7)
    segments(x0  = 1:length(sample_name),
             y0  = x$best_estimate - x$std_error,
             x1  = 1:length(sample_name),
             y1  = x$best_estimate + x$std_error,
             col = "gray")
  }
  par(old.par)
}

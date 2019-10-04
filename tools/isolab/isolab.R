#' wl-20-02-2018, Tue: commence
#' wl-24-02-2018, Sat: not install. Directly use package source code
#' wl-03-04-2018, Tue: abandon `xcms` and use directly `tsv` data formats
#' wl-10-04-2018, Tue: substantial changes
#' wl-11-04-2018, Wed: command line and change plot function
#' wl-12-04-2018, Thu: deal with empty bottom in target produced by Galaxy
#' wl-22-05-2018, Tue: test new group
#' wl-05-02-2019, Tue: load group info from file
#' wl-23-03-2019, Sat: tidy up and use Vim's folding as outline view
#' wl-24-03-2019, Sun: apply planemo test and xlsx file, prodecuded by
#'   'WriteXLS' in local is fine. So drop R package 'writexl' which cannot
#'    output row names of data frame.
#' wl-26-03-2019, Tue: drop R package 'ecipex' and use its scripts directly
#'  since no 'r-ecipex' in any conda repositories. Also extract function
#'  'makeup.R', which is called in 'ecipex' from R package 'CHNOSZ'.

## ==== General settings ====

rm(list=ls(all=T))

#' flag for command-line use or not. If false, only for debug interactively.
com_f <- T

#' galaxy will stop even if R has warning message
options(warn=-1) #' disable R warning. Turn back: options(warn=0)

#' Setup R error handling to go to stderr
options( show.error.messages=F, error = function (){
  cat( geterrmessage(), file=stderr() )
  q( "no", 1, F )
})

#' we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(optparse)
  library(WriteXLS)
  #' library(ecipex)     #' use its R scripts directly
  library(gsubfn)        #' Only for function strapply
})

## ==== Command line or interactive setting ====

if(com_f){

  #' Setup home directory
  #' wl-24-11-2017, Fri: A dummy function for the base directory. The reason
  #' to write such a function is to keep the returned values by
  #' 'commandArgs' with 'trailingOnly = FALSE' in a local environment
  #' otherwise 'parse_args' will use the results of
  #' 'commandArgs(trailingOnly = FALSE)' even with 'args =
  #' commandArgs(trailingOnly = TRUE)' in its argument area.
  func <- function(){
    argv <- commandArgs(trailingOnly = FALSE)
    path <- sub("--file=","",argv[grep("--file=",argv)])
  }
  #' prog_name <- basename(func())
  tool_dir <- paste0(dirname(func()),"/")

  option_list <-
    list(
        make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
                    help="Print extra output [default]"),
        make_option(c("-q", "--quietly"), action="store_false",
                    dest="verbose", help="Print little output"),

        #' input files
        make_option("--peak_file", type="character",
                    help="Peak table with m/z, retention time and intensity"),
        make_option("--targ_file", type="character",
                    help="Parameter data metrix in which each row is one
                          instance of setting for analysis"),

        #' group abundance estimate
        make_option("--grp", type="logical", default=TRUE,
                    help="Apply group estimates of abundance"),
        make_option("--grp_file_sel", type="character", default="yes",
                    help="Load info for group abundance estimates from file or not"),
        make_option("--grp_file", type="character",
                    help="Info file for group abundance estimates"),
        make_option("--groups", type="character", default="",
                    help="Sample group information. Delimited by commas."),

        #' plot output
        make_option("--pattern_plot", type="logical", default=TRUE,
                    help="Plot patterns"),
        make_option("--residual_plot", type="logical", default=TRUE,
                    help="Plot residuals"),
        make_option("--result_plot", type="logical", default=TRUE,
                    help="Plot results"),

        #' pdf files
        make_option("--pattern_file",type="character", default="pattern.pdf",
                    help="Save pattern plot"),
        make_option("--residual_file",type="character", default="residual.pdf",
                    help="Save residual plot"),
        make_option("--result_file",type="character", default="result.pdf",
                    help="Save result plot"),

        #' Excel files
        make_option("--summary_file",type="character",default="summary.xlsx",
                    help="Save summary results in Excel"),
        make_option("--summary_grp_file",type="character",
                    default="summary_grp.xlsx",
                    help="Save group summary results")
    )

  opt <- parse_args(object=OptionParser(option_list=option_list),
                    args = commandArgs(trailingOnly = TRUE))

} else {
  #' tool_dir <- "C:/R_lwc/isolab/"         #' for windows
  tool_dir <- "~/my_galaxy/isolab/"  #' for linux. must be case-sensitive
  opt  <- list(
      #' input files and group abundance estimate
      peak_file    = paste0(tool_dir,"test-data/xcms.tsv"),
      targ_file    = paste0(tool_dir,"test-data/xcms_tar.tsv"),
      grp          = "TRUE",
      grp_file_sel = "no",
      #' grp_file     = paste0(tool_dir,"test-data/xcms_grp.tsv"),
      groups       = "C12,C12,C12,C12,C13,C13,C13,C13",
      #' peak_file = paste0(tool_dir,"test-data/ecamam12.tsv"),
      #' targ_file = paste0(tool_dir,"test-data/ecamam12_tar.tsv"),
      #' grp_file  = paste0(tool_dir,"test-data/ecamam12_grp.tsv"),
      #' groups    = "12C_Lys,12C_Lys,12C_Lys,12C_Glu,12C_Glu,12C_Glu,12C_Lys,12C_Lys,12C_Lys,13C_Lys,13C_Lys,13C_Lys,12C_Glu,12C_Glu,12C_Glu,13C_Glu,13C_Glu,13C_Glu,12C_Lys,12C_Lys,12C_Lys,13C_Lys,13C_Lys,13C_Lys,12C_Glu,12C_Glu,12C_Glu,13C_Glu,13C_Glu,13C_Glu,12C_Lys,12C_Lys,12C_Lys,13C_Lys,13C_Lys, 13C_Lys,12C_Glu,12C_Glu,12C_Glu,13C_Glu,13C_Glu,13C_Glu",

      #' plot output
      pattern_plot  = TRUE,
      residual_plot = TRUE,
      result_plot   = TRUE,

      #' pdf files
      pattern_file  = paste0(tool_dir,"test-data/res/pattern.pdf"),
      residual_file = paste0(tool_dir,"test-data/res/residual.pdf"),
      result_file   = paste0(tool_dir,"test-data/res/result.pdf"),

      #' Excel files
      summary_file     = paste0(tool_dir,"test-data/res/summary.xlsx"),
      summary_grp_file = paste0(tool_dir,"test-data/res/summary_grp.xlsx")
  )
}

#' print(opt)

suppressPackageStartupMessages({
  source(paste0(tool_dir,"pkgs/all_IsotopicLabelling.R"))
  source(paste0(tool_dir,"pkgs/all_ecipex.R"))
  load(paste0(tool_dir,"pkgs/nistiso.rda"))  #' needed by ecipex
})


## ==== 1) Data preparation ====

#' Load peak table
peak <- read.table(opt$peak_file, header = T, sep = "\t",
                   fill = T,stringsAsFactors = F)
#' cat("\n the row names is\n")

#' Load batch parameters
targets <- read.table(opt$targ_file, header = T, sep = "\t",
                      fill = T,stringsAsFactors = F)

#' wl-13-04-2018, Fri: targets file produced by galaxy will have empty
#' bottom rows. Must remove any empty rows.
tmp     <- targets[,-ncol(targets)]
idx     <- complete.cases(tmp)
targets <- targets[idx,]

#' transpose data frame
targets  <- as.data.frame(t(targets))
#' targets  <- as.data.frame(t(targets),stringsAsFactors = F)
#' targets  <- as.list(targets)

## ==== 2) Batch process ====

res_bat <- lapply(targets,function(x){  #'  x = targets[[1]]

  info <- isotopic_information(compound  = as.character(x["compound"]),
                               charge    = as.numeric(as.character(x["charge"])),
                               labelling = as.character(x["labelling"]))

  patterns <- isotopic_pattern(peak_table  = peak,
                               info        = info,
                               mass_shift  = as.numeric(as.character(x["mass_shift"])),
                               RT          = as.numeric(as.character(x["RT"])),
                               RT_shift    = as.numeric(as.character(x["RT_shift"])),
                               chrom_width = as.numeric(as.character(x["chrom_width"])))

  res <- find_abundance(patterns          = patterns,
                        info              = info,
                        initial_abundance = as.numeric(as.character(x["initial_abundance"])),
                        charge            = as.numeric(as.character(x["charge"])))

  return(res)
})
names(res_bat) <- as.character(unlist(targets[2,]))

## ==== 3) Process the results ====

#' Plots
if (opt$pattern_plot) {
  pdf(file = opt$pattern_file, onefile = T,width=15, height=10)
  lapply(res_bat, function(x) plot.func(x=x, type="patterns"))
  dev.off()
}

if (opt$residual_plot) {
  pdf(file = opt$residual_file, onefile = T,width=15, height=10)
  lapply(res_bat, function(x) plot.func(x=x, type="residuals"))
  dev.off()
}

if (opt$result_plot) {
  pdf(file = opt$result_file, onefile = T,width=15, height=10)
  lapply(res_bat, function(x) plot.func(x=x, type="summary"))
  dev.off()
}

#' Summary of individual sample
summ  <- lapply(res_bat, function(x) as.data.frame(summary(x)))
WriteXLS(summ, ExcelFileName = opt$summary_file, row.names = T, FreezeRow = 1)

#' Summary of grouped samples
if (opt$grp) {
  #' get group info
  if (opt$grp_file_sel == "yes") {
    groups <- read.table(opt$grp_file, header = FALSE, sep = "\t",
                         stringsAsFactors = F)
    groups <- groups[,1,drop = TRUE]
    #' wl-30-11-2018, Fri: group file must be one column without header. The 
    #'  file extension can be tsv, csv or txt. sep="\t" takes no effect on one 
    #'  column file.
  } else {
    groups <- opt$groups
    groups <- unlist(strsplit(groups,","))
    groups <- gsub("^[ \t]+|[ \t]+$", "", groups)  #' trim white spaces
  }
  groups <- as.factor(tolower(groups))

  summ_grp <- lapply(res_bat, function(x) group_labelling(x, groups = groups))
  WriteXLS(summ_grp, ExcelFileName = opt$summary_grp_file, row.names = T, FreezeRow = 1)
}

## ==== DEBUG: Step-by-step codes ====

#' wl-13-04-2018, Fri: test our own data set

if (F) {
  info <- isotopic_information(compound="X51H98O6", labelling="C")
  #' names(info)
  #' info$isotopes

  patterns <- isotopic_pattern(peak, info, mass_shift=0.05,
                               RT=385, RT_shift=10, chrom_width=9)
  #' View(patterns)

  fitted <- find_abundance(patterns=patterns, info=info,
                           initial_abundance=NA, charge=1)
  #' names(fitted)
  #' summary(fitted)
  #' save_labelling(fitted)

  #' Or use wrapper function
  fitted <- main_labelling(peak, compound="X51H98O6",
                           charge=1, labelling="C", mass_shift=0.05,
                           RT=380, RT_shift=10, chrom_width=7,
                           initial_abundance=NA)

  plot(x=fitted, type="patterns", saveplots=F)
  plot(x=fitted, type="residuals", saveplots=F)
  plot(x=fitted, type="summary", saveplots=F)

  #' Group the samples and obtain grouped estimates
  #' groups <- factor(c(rep("C12",4), rep("C13",4))) 
  groups  <- "12C_Lys,12C_Lys,12C_Lys,12C_Glu,12C_Glu,12C_Glu,12C_Lys, 
              12C_Lys,12C_Lys,13C_Lys,13C_Lys,13C_Lys,12C_Glu,12C_Glu, 
              12C_Glu,13C_Glu,13C_Glu,13C_Glu,12C_Lys,12C_Lys,12C_Lys, 
              13C_Lys,13C_Lys,13C_Lys,12C_Glu,12C_Glu,12C_Glu,13C_Glu, 
              13C_Glu,13C_Glu,12C_Lys,12C_Lys,12C_Lys,13C_Lys,13C_Lys, 
              13C_Lys,12C_Glu,12C_Glu,12C_Glu,13C_Glu,13C_Glu,13C_Glu"
  groups  <- unlist(strsplit(groups,","))
  groups  <- gsub("^[ \t]+|[ \t]+$", "", groups)  #' trim white spaces
  groups  <- factor(groups)

  grp_est <- group_labelling(fitted,groups=groups)
  grp_est
}

## ==== DEBUG: Original batch process ====

if (F) {
  #' Batch-process
  #' wl-13-04-2018, Fri: Not work. Here 'targets' should not be transposed.
  #' wl-30-05-2018, Wed: Should change the original R codes

  bat_grp_est <- batch_labelling(peak_table=peak, targets=targets,
                                 groups=groups,
                                 plot_patterns=F, plot_residuals=F,
                                 plot_results=F, save_results=F)
  bat_grp_est

  #' Use data set from xcms
  load("./test-data/xcms_obj.rda") #' data("xcms_obj")
  peak <- table_xcms(xcms_obj)
  write.table(peak, file="./test-data/xcms_obj.tsv", sep = "\t",
              row.names = FALSE, quote = FALSE)

  #' Get the example data frame containing target abalytes
  load("./test-data/targets.rda") #' data("targets")
  write.table(targets,file="./test-data/targets.tsv", sep="\t",
              row.names = FALSE, quote = FALSE)

  para  <- c("compound", "charge", "labelling", "RT", "RT_shift",
             "chrom_width", "mass_shift", "initial_abundance")
  targets <- targets[,para]

  res <- main_labelling(peak_table        = peak,
                        compound          = as.character(x["compound"]),
                        charge            = as.numeric(as.character(x["charge"])),
                        labelling         = as.character(x["labelling"]),
                        mass_shift        = as.numeric(as.character(x["mass_shift"])),
                        RT                = as.numeric(as.character(x["RT"])),
                        RT_shift          = as.numeric(as.character(x["RT_shift"])),
                        chrom_width       = as.numeric(as.character(x["chrom_width"])),
                        initial_abundance = as.numeric(as.character(x["initial_abundance"])))

  summ     <- lapply(res_bat,"[[","summary")
}

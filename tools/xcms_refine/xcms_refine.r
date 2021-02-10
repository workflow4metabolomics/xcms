#!/usr/bin/env Rscript

# ----- LOG FILE -----
log_file=file("log.txt", open = "wt")
sink(log_file)
sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname){ argv <- commandArgs(trailingOnly=FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep="/")) }
source_local("lib.r")

pkgs <- c("xcms","batch","RColorBrewer")
loadAndDisplayPackages(pkgs)
cat("\n\n");

# ----- EXPLORE NAMESPACE -----
# cat("\tXCMS NAMESPACE INFO\n")
# ls(asNamespace("xcms"))

# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args <- parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names=F, quote=F, sep='\t')

cat("\n\n")

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

#saving the specific parameters
args_method  <- args$method
args_image   <- args$image
args_msLevel <- args$msLevel
#args[c('method','image','msLevel')] <- NULL
param_args <- list()

if ( args_method == "CleanPeaks" ) {
  param_args$maxPeakwidth <- as.numeric(args$maxPeakwidth)
  if (is.na(as.numeric(param_args$maxPeakwidth))) stop("\n\nERROR: The maxPeakwidth argument cannot be coerced to a numeric value.")
} else if ( args_method == "FilterIntensity" ) {
  param_args$threshold <- as.numeric(args$threshold)
  if (is.na(as.numeric(param_args$threshold))) stop("\n\nERROR: The threshold argument cannot be coerced to a numeric value.")
  param_args$nValues <- as.numeric(args$nValues)
  if (is.na(as.numeric(param_args$nValues))) stop("\n\nERROR: The nValues argument cannot be coerced to a numeric value.")
  if (as.integer(param_args$nValues) != param_args$nValues) stop("\n\nERROR: The nValues argument is not an integer value.")
  param_args$value <- args$value
} else if ( args_method == "MergeNeighboringPeaks" ) {
  if (is.na(as.numeric(args$expandRt))) stop("\n\nERROR: The expandRt argument cannot be coerced to a numeric value.")
  if (is.na(as.numeric(args$expandMz))) stop("\n\nERROR: The expandMz argument cannot be coerced to a numeric value.")
  if (is.na(as.numeric(args$ppm))) stop("\n\nERROR: The ppm argument cannot be coerced to a numeric value.")
  if (is.na(as.numeric(args$minProp))) stop("\n\nERROR: The minProp argument cannot be coerced to a numeric value.")
  param_args$expandRt <- args$expandRt
  param_args$expandMz <- args$expandMz
  param_args$ppm      <- args$ppm
  param_args$minProp  <- args$minProp
}

cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

#image is an .RData file necessary to use xset variable given by previous tools
load(args_image)
if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. Such RData as this might have been created by an old version of XMCS 2.*")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
rawFilePath <- retrieveRawfileInTheWorkingDirectory(singlefile, zipfile, args)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile

cat("\n\n")


# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


cat("\t\tPREPARE PARAMETERS\n\n")

if ( args_method == "CleanPeaks" ) {
  refineChromPeaksParam <- CleanPeaksParam(maxPeakwidth = param_args$maxPeakwidth)
} else if ( args_method == "FilterIntensity" ) {
  refineChromPeaksParam <- FilterIntensityParam(
    threshold = param_args$threshold,
    nValues = param_args$nValues,
    value = param_args$value
  )
} else if ( args_method == "MergeNeighboringPeaks" ) {
  refineChromPeaksParam <- MergeNeighboringPeaksParam(
    expandRt = param_args$expandRt,
    expandMz = param_args$expandMz,
    ppm = param_args$ppm,
    minProp = param_args$minProp
  )
}

cat(str(refineChromPeaksParam))

cat("\n\n\t\tCOMPUTE\n")

# cat("\t\t\tACE DEBUG show the xdata\n")
# cat(str(xdata))
# cat("\n\t\t\tACE DEBUG (end show the xdata)\n")

# cat("\t\t\tACE DEBUG show updateObject(xdata)\n")
xdata <- updateObject(xdata)
# cat(str(xdata))
# cat("\n\t\t\tACE DEBUG (end show updateObject(xdata))\n")
# 
# cat("\t\t\tACE DEBUG show the param_args\n")
# cat(str(param_args))
# cat("\n\t\t\tACE DEBUG (end show the param_args)\n")
# 
# cat("\t\t\tACE DEBUG show the refineChromPeaksParam\n")
# cat(str(refineChromPeaksParam))
# cat("\n\t\t\tACE DEBUG (end show the refineChromPeaksParam)\n")
# 
# cat("\t\t\tACE DEBUG show the slotNames\n")
# cat(str(slotNames(refineChromPeaksParam)))
# cat("\n\t\t\tACE DEBUG (end show the slotNames)\n")

# # clear the arguement list to remove unexpected key/value as singlefile_galaxyPath or args_method ...
# param_args <- param_args[names(param_args) %in% slotNames(do.call(paste0(args_method,"Param"), param_args))]
#
# refineChromPeaksParam <- do.call(paste0(args_method,"Param"), param_args)
# print(refineChromPeaksParam)

# cat("\t\t\tACE DEBUG invoke refineChromPeaks\n")
xdata <- refineChromPeaks(xdata, param = refineChromPeaksParam)
# cat("\t\t\tACE DEBUG returned from refineChromPeaks\n")
#
#
# cat("\t\tDRAW GRAPHICS\n")
# getPlotChromPeakDensity(xdata, param = groupChromPeaksParam)
#
# if (exists("intval")) {
#     getPeaklistW4M(xdata, intval, convertRTMinute, numDigitsMZ, numDigitsRT, naTOzero, "variableMetadata.tsv", "dataMatrix.tsv")
# }
#
cat("\n\n")

# ----- EXPORT -----

cat("\tXCMSnExp OBJECT INFO\n")
print(xdata)
cat("\n\n")

cat("\txcmsSet OBJECT INFO\n")
# Get the legacy xcmsSet object
xset <- getxcmsSetObject(xdata)
print(xset)
cat("\n\n")

#saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xdata", "zipfile", "singlefile", "md5sumList", "sampleNamesList")
save(list=objects2save[objects2save %in% ls()], file="xcmsSet.RData")

cat("\n\n")


cat("\tDONE\n")

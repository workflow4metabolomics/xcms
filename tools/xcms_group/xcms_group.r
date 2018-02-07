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

pkgs=c("xcms","batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
listArguments <- parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(listArguments), col.names=F, quote=F, sep='\t')

cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

#image is an .RData file necessary to use xset variable given by previous tools
load(listArguments[["image"]]); listArguments[["image"]]=NULL
if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. This RData should have been created by an old version of XMCS 2.*")

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, listArguments)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
listArguments <- rawFilePath$listArguments
directory <- retrieveRawfileInTheWorkingDirectory(singlefile, zipfile)

# Check some character issues
md5sumList <- list("origin" = getMd5sum(directory))
checkXmlStructure(directory)
checkFilesCompatibilityWithXcms(directory)


cat("\n\n")

#Import the different functions

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

# Save arguments to generate a report
if (!exists("listOFlistArguments")) listOFlistArguments <- list()
listOFlistArguments[[paste(format(Sys.time(), "%y%m%d-%H:%M:%S_"),listArguments[["xfunction"]],sep="")]] <- listArguments


#saving the commun parameters
BPPARAM <- MulticoreParam(1)
if (!is.null(listArguments[["BPPARAM"]])){
    BPPARAM <- MulticoreParam(listArguments[["BPPARAM"]]); listArguments[["BPPARAM"]] <- NULL
}
register(BPPARAM)

#saving the specific parameters
method <- listArguments[["method"]]; listArguments[["method"]] <- NULL

if (!is.null(listArguments[["convertRTMinute"]])){
    convertRTMinute <- listArguments[["convertRTMinute"]]; listArguments[["convertRTMinute"]] <- NULL
}
if (!is.null(listArguments[["numDigitsMZ"]])){
    numDigitsMZ <- listArguments[["numDigitsMZ"]]; listArguments[["numDigitsMZ"]] <- NULL
}
if (!is.null(listArguments[["numDigitsRT"]])){
    numDigitsRT <- listArguments[["numDigitsRT"]]; listArguments[["numDigitsRT"]] <- NULL
}
if (!is.null(listArguments[["intval"]])){
    intval <- listArguments[["intval"]]; listArguments[["intval"]] <- NULL
}

cat("\n\n")


# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")

#change the default display settings
pdf(file="Rplots.pdf", width=16, height=12)
par(mfrow=c(2,2))


cat("\t\tCOMPUTE\n")

cat("\t\t\tPerform the correspondence\n")
listArguments[["sampleGroups"]] = xdata$sample_group
groupChromPeaksParam <- do.call(paste0(method,"Param"), listArguments)
print(groupChromPeaksParam)
xdata <- groupChromPeaks(xdata, param = groupChromPeaksParam)

# Get the legacy xcmsSet object
suppressWarnings(xset <- as(xdata, 'xcmsSet'))
sampclass(xset) <- xset@phenoData$sample_group

cat("\n\n")

dev.off()

getPeaklistW4M(xdata, intval, convertRTMinute, numDigitsMZ, numDigitsRT, "variableMetadata.tsv", "dataMatrix.tsv")

cat("\n\n")

# ----- EXPORT -----

cat("\tXCMSnExp OBJECT INFO\n")
print(xdata)
cat("\n\n")

cat("\txcmsSet OBJECT INFO\n")
print(xset)
cat("\n\n")

#saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xdata", "zipfile", "singlefile", "listOFlistArguments", "md5sumList", "sampleNamesList")
save(list=objects2save[objects2save %in% ls()], file="group.RData")

cat("\n\n")


cat("\tDONE\n")

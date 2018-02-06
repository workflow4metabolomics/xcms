#!/usr/bin/env Rscript

# ----- LOG FILE -----
log_file <- file("log.txt", open="wt")
sink(log_file)
sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

pkgs=c("xcms","batch")
for(pkg in pkgs) suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))

sessioninfo = sessionInfo()
cat(sessioninfo$R.version$version.string,"\n")
cat("Main packages:\n")
for (pkg in names(sessioninfo$otherPkgs)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
cat("Other loaded packages:\n")
for (pkg in names(sessioninfo$loadedOnly)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")

source_local <- function(fname){ argv <- commandArgs(trailingOnly=FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep="/")) }
cat("\n\n");


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
listArguments <- parseCommandArgs(evaluate = FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(listArguments), col.names=F, quote=F, sep='\t')

cat("\n\n");


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

#Import the different functions
source_local("lib.r")

cat("\n\n")

#Import the different functions

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

# Save arguments to generate a report
if (!exists("listOFlistArguments")) listOFlistArguments <- list()
listOFlistArguments[[paste(format(Sys.time(), "%y%m%d-%H:%M:%S_"), listArguments[["xfunction"]], sep="")]] = listArguments


#saving the commun parameters


BPPARAM = MulticoreParam(1)
if (!is.null(listArguments[["nSlaves"]])){
    BPPARAM <- MulticoreParam(listArguments[["nSlaves"]]); listArguments[["nSlaves"]] <- NULL
}
register(BPPARAM)

thefunction <- listArguments[["xfunction"]]; listArguments[["xfunction"]] <- NULL #delete from the list of arguments



#saving the specific parameters
if (!is.null(listArguments[["method"]])){
    method <- listArguments[["method"]]; listArguments[["method"]] <- NULL
}

ticspdf <- "TICs.pdf"
bicspdf <- "BICs.pdf"

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




# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


cat("\t\tCOMPUTE\n")
files <- getMSFiles(directory)

# PhenoData
s_groups <- sapply(files, function(x) tail(unlist(strsplit(dirname(x),"/")), n=1))
s_name <- file_path_sans_ext(basename(files))
pd <- data.frame(sample_name=s_name, sample_group=s_groups, stringsAsFactors=FALSE)

# Reading Raw
raw_data <- readMSData(files=files, pdata = new("NAnnotatedDataFrame", pd), mode="onDisk")

# Peak Calling
if (method == "centWave") methodParam <- "CentWaveParam"
findChromPeaksParam <- do.call(methodParam, listArguments)
xdata <- findChromPeaks(raw_data, param=findChromPeaksParam)

# check if there are no peaks
if (nrow(chromPeaks(xdata)) == 0) {
    stop("No peaks were detected. You should review your settings")
}

# transform the files absolute pathways into relative pathways
xdata@processingData@files <- sub(paste(getwd(), "/", sep="") , "", xdata@processingData@files)

# create a sampleMetada file
sampleNamesList <- getSampleMetadata(xdata=xdata, sampleMetadataOutput=sampleMetadataOutput)

cat("\n\n")


# -- TIC --
cat("\t\tGET TIC GRAPH\n")
#@TODO: one day, use xdata instead of xset to draw the TICs and BPC or a complete other method
xset <- as(xdata, 'xcmsSet')
sampclass(xset) <- xset@phenoData$sample_group

getTICs(xdata=xdata, pdfname=ticspdf, rt="raw")
getBPCs(xdata=xdata, rt="raw", pdfname=bicspdf)

cat("\n\n")

# ----- EXPORT -----

cat("\tXSET OBJECT INFO\n")
print(xdata)

#saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xdata", "zipfile", "singlefile", "listOFlistArguments", "md5sumList", "sampleNamesList")
xsetRdataOutput <- paste(thefunction, "RData", sep=".")
save(list=objects2save[objects2save %in% ls()], file=xsetRdataOutput)

cat("\n\n")


cat("\tDONE\n")

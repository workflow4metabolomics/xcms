#!/usr/bin/env Rscript

#Import the different functions
source_local <- function(fname){ argv <- commandArgs(trailingOnly=FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep="/")) }
source_local("lib.r")

pkgs <- c("xcms","batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");

listArguments <- parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects

# Handle infiles
if (!exists("singlefile")) singlefile <- NULL
if (!exists("zipfile")) zipfile <- NULL
rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, listArguments)
zipfile <- rawFilePath$zipfile
singlefile <- rawFilePath$singlefile
listArguments <- rawFilePath$args
directory <- retrieveRawfileInTheWorkingDirectory(singlefile, zipfile)

cat("\tXSET MERGING...\n")

for(image in listArguments[["images"]]) {
    load(image)
    cat(sampleNamesList$sampleNamesOrigin,"\n")
    if (!exists("xdata_merged")) {
        xdata_merged <- xdata
        singlefile_merged <- singlefile
        md5sumList_merged <- md5sumList
        sampleNamesList_merged <- sampleNamesList
    } else {
        xdata_merged <- c(xdata_merged,xdata)
        singlefile_merged <- c(singlefile_merged,singlefile)
        md5sumList_merged$origin <- rbind(md5sumList_merged$origin,md5sumList$origin)
        sampleNamesList_merged$sampleNamesOrigin <- c(sampleNamesList_merged$sampleNamesOrigin,sampleNamesList$sampleNamesOrigin)
        sampleNamesList_merged$sampleNamesMakeNames <- c(sampleNamesList_merged$sampleNamesMakeNames,sampleNamesList$sampleNamesMakeNames)
    }
}
rm(image)
xdata <- xdata_merged; rm(xdata_merged)
singlefile <- singlefile_merged; rm(singlefile_merged)
md5sumList <- md5sumList_merged; rm(md5sumList_merged)
sampleNamesList <- sampleNamesList_merged; rm(sampleNamesList_merged)

if (!is.null(listArguments[["sampleMetadata"]])) {
    cat("\tXSET PHENODATA SETTING...\n")
    sampleMetadataFile <- listArguments[["sampleMetadata"]]
    sampleMetadata <- read.table(sampleMetadataFile, h=F, sep=";", stringsAsFactors=F)
    if (ncol(sampleMetadata) < 2) sampleMetadata <- read.table(sampleMetadataFile, h=F, sep="\t", stringsAsFactors=F)
    if (ncol(sampleMetadata) < 2) sampleMetadata <- read.table(sampleMetadataFile, h=F, sep=",", stringsAsFactors=F)
    if (ncol(sampleMetadata) < 2) {
        error_message="Your sampleMetadata file seems not well formatted. The column separators accepted are ; , and tabulation"
        print(error_message)
        stop(error_message)
    }
    xdata@phenoData@data$sample_group=sampleMetadata$V2[match(xdata@phenoData@data$sample_name,sampleMetadata$V1)]

    if (any(is.na(pData(xdata)$sample_group))) {
        sample_missing <- pData(xdata)$sample_name[is.na(pData(xdata)$sample_group)]
        error_message <- paste("Those samples are missing in your sampleMetadata:", paste(sample_missing, collapse=" "))
        print(error_message)
        stop(error_message)
    }
}
save.image()

cat("\tXCMSnExp OBJECT INFO\n")
print(pData(xdata))
print(xdata)
cat("\n\n")

cat("\txcmsSet OBJECT INFO\n")
# Get the legacy xcmsSet object
xset <- getxcmsSetObject(xdata)
print(xset@phenoData)
print(xset)
cat("\n\n")

cat("\tSAVE RData\n")
#saving R data in .Rdata file to save the variables used in the present tool
objects2save <- c("xdata", "zipfile", "singlefile", "md5sumList", "sampleNamesList")
save(list=objects2save[objects2save %in% ls()], file="merged.RData")

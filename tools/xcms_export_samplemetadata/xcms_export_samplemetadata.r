#!/usr/bin/env Rscript

#Import the different functions
source_local <- function(fname){ argv <- commandArgs(trailingOnly=FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep="/")) }
source_local("lib.r")
source_local("lib-xcms3.x.x.r")

pkgs <- c("xcms","batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");

args <- parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects


sampleMetadata <- NULL
for(image in args$images) {
    load(image)
    if (is.null(sampleMetadata))
        sampleMetadata <- xdata@phenoData@data
    else
        sampleMetadata <- rbind(sampleMetadata,xdata@phenoData@data)
}
colnames(sampleMetadata) <- c("sample_name","class")
sampleMetadata$sample_name <- make.names(sampleNamesOrigin)

# Create a sampleMetada file
write.table(sampleMetadata,file="sampleMetadata.tsv", sep="\t", row.names=FALSE, quote=FALSE)

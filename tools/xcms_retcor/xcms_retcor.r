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

pkgs <- c("xcms","batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
listArguments = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(listArguments), col.names=F, quote=F, sep='\t')

cat("\n\n")

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

#saving the commun parameters
BPPARAM = MulticoreParam(1)
if (!is.null(listArguments[["nSlaves"]])){
    BPPARAM = MulticoreParam(listArguments[["nSlaves"]]); listArguments[["nSlaves"]]=NULL
}
register(BPPARAM)

#saving the specific parameters
method <- listArguments[["method"]]; listArguments[["method"]] <- NULL

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
listArguments <- rawFilePath$args
directory <- retrieveRawfileInTheWorkingDirectory(singlefile, zipfile)

# Check some character issues
md5sumList <- list("origin" = getMd5sum(directory))
checkXmlStructure(directory)
checkFilesCompatibilityWithXcms(directory)


cat("\n\n")


# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


cat("\t\tCOMPUTE\n")

#change the default display settings
pdf(file="Rplots.pdf", width=16, height=12)
#try to change the legend display
#par(xpd=NA)
#par(xpd=T, mar=par()$mar+c(0,0,0,4))

cat("\t\t\tAlignment/Retention Time correction\n")
adjustRtimeParam <- do.call(paste0(method,"Param"), listArguments)
print(adjustRtimeParam)
xdata <- adjustRtime(xdata, param=adjustRtimeParam)

dev.off() #dev.new(file="Rplots.pdf", width=16, height=12)# Get the legacy xcmsSet object

# Get the legacy xcmsSet object
xset <- getxcmsSetObject(xdata)

cat("\n\n")


# -- TIC --
cat("\t\tGET TIC GRAPH\n")
#@TODO: one day, use xdata instead of xset to draw the TICs and BPC or a complete other method
getTICs(xcmsSet=xset, rt="raw", pdfname="TICs.pdf")
getBPCs(xcmsSet=xset, rt="raw", pdfname="BICs.pdf")

cat("\n\n")

# ----- EXPORT -----

cat("\tXCMSnExp OBJECT INFO\n")
print(xdata)
cat("\n\n")

cat("\txcmsSet OBJECT INFO\n")
print(xset)
cat("\n\n")

#saving R data in .Rdata file to save the variables used in the present tool
objects2save = c("xdata","zipfile","singlefile","md5sumList","sampleNamesList")
save(list=objects2save[objects2save %in% ls()], file="retcor.RData")

cat("\n\n")


cat("\tDONE\n")

#!/usr/bin/env Rscript
# CAMERA.r version="2.2.1"



# ----- PACKAGE -----
cat("\tPACKAGE INFO\n")

setRepositories(graphics=F, ind=31)

#pkgs=c("xcms","batch")
pkgs=c("parallel","BiocGenerics", "Biobase", "Rcpp", "mzR", "xcms","snow","igraph","CAMERA","multtest","batch")
for(p in pkgs) {
    suppressPackageStartupMessages(suppressWarnings(library(p, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))
    cat(p,"\t",as.character(packageVersion(p)),"\n",sep="")
}
source_local <- function(fname){ argv <- commandArgs(trailingOnly = FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep="/")) }

cat("\n\n");



# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")

listArguments = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(listArguments), col.names=F, quote=F, sep='\t')

cat("\n\n");


# ----- PROCESSING INFILE -----
cat("\tINFILE PROCESSING INFO\n")

#image is an .RData file necessary to use xset variable given by previous tools
if (!is.null(listArguments[["image"]])){
    load(listArguments[["image"]]); listArguments[["image"]]=NULL
}

if (listArguments[["xfunction"]] %in% c("combinexsAnnos")) {
    load(listArguments[["image_pos"]])
    xaP=xa
    listOFlistArgumentsP=listOFlistArguments
    if (exists("xsAnnotate_object")) xaP=xsAnnotate_object

    diffrepP=NULL
    if (exists("diffrep")) diffrepP=diffrep

    load(listArguments[["image_neg"]])
    xaN=xa
    listOFlistArgumentsN=listOFlistArguments
    if (exists("xsAnnotate_object")) xaN=xsAnnotate_object

    diffrepN=NULL
    if (exists("diffrep")) diffrepN=diffrep
}


cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tARGUMENTS PROCESSING INFO\n")

# Save arguments to generate a report
if (!exists("listOFlistArguments")) listOFlistArguments=list()
listOFlistArguments[[paste(format(Sys.time(), "%y%m%d-%H:%M:%S_"),listArguments[["xfunction"]],sep="")]] = listArguments


#saving the commun parameters
thefunction = listArguments[["xfunction"]]
listArguments[["xfunction"]]=NULL #delete from the list of arguments

xsetRdataOutput = paste(thefunction,"RData",sep=".")
if (!is.null(listArguments[["xsetRdataOutput"]])){
    xsetRdataOutput = listArguments[["xsetRdataOutput"]]; listArguments[["xsetRdataOutput"]]=NULL
}

rplotspdf = "Rplots.pdf"
if (!is.null(listArguments[["rplotspdf"]])){
    rplotspdf = listArguments[["rplotspdf"]]; listArguments[["rplotspdf"]]=NULL
}

dataMatrixOutput = "dataMatrix.tsv"
if (!is.null(listArguments[["dataMatrixOutput"]])){
    dataMatrixOutput = listArguments[["dataMatrixOutput"]]; listArguments[["dataMatrixOutput"]]=NULL
}

variableMetadataOutput = "variableMetadata.tsv"
if (!is.null(listArguments[["variableMetadataOutput"]])){
    variableMetadataOutput = listArguments[["variableMetadataOutput"]]; listArguments[["variableMetadataOutput"]]=NULL
}

if (!is.null(listArguments[["new_file_path"]])){
    new_file_path = listArguments[["new_file_path"]]; listArguments[["new_file_path"]]=NULL
}

#Import the different functions
source_local("lib.r")

# We unzip automatically the chromatograms from the zip files.
if (thefunction %in% c("annotatediff"))  {
    rawFilePath = getRawfilePathFromArguments(listArguments)
    zipfile = rawFilePath$zipfile
    singlefile = rawFilePath$singlefile
    listArguments = rawFilePath$listArguments
    directory = retrieveRawfileInTheWorkingDirectory(singlefile, zipfile)
}


#addition of xset object to the list of arguments in the first position
if (exists("xset") != 0){
    listArguments=append(list(xset), listArguments)
}

cat("\n\n")




# ----- PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")

#change the default display settings
pdf(file=rplotspdf, width=16, height=12)


if (thefunction %in% c("annotatediff")) {
    results_list=annotatediff(xset=xset,listArguments=listArguments,variableMetadataOutput=variableMetadataOutput,dataMatrixOutput=dataMatrixOutput,new_file_path=new_file_path)
    xa=results_list[["xa"]]
    diffrep=results_list[["diffrep"]]
    variableMetadata=results_list[["variableMetadata"]]

    cat("\n\n")
    cat("\tXSET OBJECT INFO\n")
    print(xa)
}

if (thefunction %in% c("combinexsAnnos")) {
    cAnnot=combinexsAnnos_function(
        xaP=xaP,xaN=xaN,
        listOFlistArgumentsP=listOFlistArgumentsP,listOFlistArgumentsN=listOFlistArgumentsN,
        diffrepP=diffrepP,diffrepN=diffrepN,
        pos=listArguments[["pos"]],tol=listArguments[["tol"]],ruleset=listArguments[["ruleset"]],keep_meta=listArguments[["keep_meta"]],
        intval=listArguments[["intval"]], convertRTMinute=listArguments[["convertRTMinute"]], numDigitsMZ=listArguments[["numDigitsMZ"]], numDigitsRT=listArguments[["numDigitsRT"]],
        variableMetadataOutput=variableMetadataOutput
    )
}

dev.off()


#saving R data in .Rdata file to save the variables used in the present tool
objects2save = c("xa","variableMetadata","diffrep","cAnnot","listOFlistArguments","zipfile","singlefile")
save(list=objects2save[objects2save %in% ls()], file=xsetRdataOutput)

cat("\n\n")

cat("\tDONE\n")

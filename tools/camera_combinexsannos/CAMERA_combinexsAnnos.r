#!/usr/bin/env Rscript
# CAMERA.r version="2.2.1"



# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname){ argv <- commandArgs(trailingOnly=FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep="/")) }
source_local("lib.r")

pkgs=c("CAMERA","multtest","batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");



# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")

args = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names=F, quote=F, sep='\t')

cat("\n\n");


# ----- PROCESSING INFILE -----
cat("\tINFILE PROCESSING INFO\n")

#image is an .RData file necessary to use xset variable given by previous tools
if (!is.null(args$image)){
    load(args$image); args$image=NULL
}

if (args$xfunction %in% c("combinexsAnnos")) {
    load(args$image_pos)
    xaP=xa
    if (exists("xsAnnotate_object")) xaP=xsAnnotate_object

    diffrepP=NULL
    if (exists("diffrep")) diffrepP=diffrep

    load(args$image_neg)
    xaN=xa
    if (exists("xsAnnotate_object")) xaN=xsAnnotate_object

    diffrepN=NULL
    if (exists("diffrep")) diffrepN=diffrep
}


cat("\n\n")


# ----- ARGUMENTS PROCESSING -----
cat("\tARGUMENTS PROCESSING INFO\n")

# Save arguments to generate a report
if (!exists("listOFargs")) listOFargs=list()
listOFargs[[paste(format(Sys.time(), "%y%m%d-%H:%M:%S_"),args$xfunction,sep="")]] = args


#saving the commun parameters
thefunction = args$xfunction
args$xfunction=NULL #delete from the list of arguments

xsetRdataOutput = paste(thefunction,"RData",sep=".")
if (!is.null(args$xsetRdataOutput)){
    xsetRdataOutput = args$xsetRdataOutput; args$xsetRdataOutput=NULL
}

rplotspdf = "Rplots.pdf"
if (!is.null(args$rplotspdf)){
    rplotspdf = args$rplotspdf; args$rplotspdf=NULL
}

variableMetadataOutput = "variableMetadata.tsv"
if (!is.null(args$variableMetadataOutput)){
    variableMetadataOutput = args$variableMetadataOutput; args$variableMetadataOutput=NULL
}

# We unzip automatically the chromatograms from the zip files.
if (thefunction %in% c("annotatediff"))  {
    if (!exists("zipfile")) zipfile=NULL
    if (!exists("singlefile")) singlefile=NULL
    rawFilePath = getRawfilePathFromArguments(singlefile, zipfile, args)
    zipfile = rawFilePath$zipfile
    singlefile = rawFilePath$singlefile
    args = rawFilePath$args
    directory = retrieveRawfileInTheWorkingDirectory(singlefile, zipfile)
}

# Because so far CAMERA isn't compatible with the new XCMSnExp object
if (exists("xdata")){
    xset <- getxcmsSetObject(xdata)
}

# addition of xset object to the list of arguments in the first position
if (exists("xset")){
    args=append(list(xset), args)
}

cat("\n\n")




# ----- PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")

#change the default display settings
pdf(file=rplotspdf, width=16, height=12)

if (thefunction %in% c("annotatediff")) {
    results_list=annotatediff(xset=xset,args=args,variableMetadataOutput=variableMetadataOutput)
    xa=results_list$xa
    diffrep=results_list$diffrep
    variableMetadata=results_list$variableMetadata

    cat("\n\n")
    cat("\tXSET OBJECT INFO\n")
    print(xa)
}

if (thefunction %in% c("combinexsAnnos")) {
    cAnnot=combinexsAnnos_function(
        xaP=xaP,xaN=xaN,
        diffrepP=diffrepP, diffrepN=diffrepN,
        pos=args$pos, tol=args$tol,ruleset=args$ruleset, keep_meta=args$keep_meta,
        convertRTMinute=args$convertRTMinute, numDigitsMZ=args$numDigitsMZ, numDigitsRT=args$numDigitsRT,
        variableMetadataOutput=variableMetadataOutput
    )
}

dev.off()


#saving R data in .Rdata file to save the variables used in the present tool
objects2save = c("xa","variableMetadata","diffrep","cAnnot","listOFargs","zipfile","singlefile")
save(list=objects2save[objects2save %in% ls()], file=xsetRdataOutput)

cat("\n\n")

cat("\tDONE\n")

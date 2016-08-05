#!/usr/bin/env Rscript
# xcms.r version="2.2.0"
#Authors ABIMS TEAM
#BPC Addition from Y.guitton


# ----- LOG FILE -----
log_file=file("log.txt", open = "wt")
sink(log_file)
sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tPACKAGE INFO\n")
#pkgs=c("xcms","batch")
pkgs=c("parallel","BiocGenerics", "Biobase", "Rcpp", "mzR", "xcms","snow","batch")
for(pkg in pkgs) {
  suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))
  cat(pkg,"\t",as.character(packageVersion(pkg)),"\n",sep="")
}
source_local <- function(fname){ argv <- commandArgs(trailingOnly = FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep="/")) }
cat("\n\n");





# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
listArguments = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(listArguments), col.names=F, quote=F, sep='\t')

cat("\n\n");


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")

#image is an .RData file necessary to use xset variable given by previous tools
if (!is.null(listArguments[["image"]])){
  load(listArguments[["image"]]); listArguments[["image"]]=NULL
}

#Import the different functions
source_local("lib.r")

cat("\n\n")

#Import the different functions

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")

# Save arguments to generate a report
if (!exists("listOFlistArguments")) listOFlistArguments=list()
listOFlistArguments[[paste(format(Sys.time(), "%y%m%d-%H:%M:%S_"),listArguments[["xfunction"]],sep="")]] = listArguments


#saving the commun parameters
thefunction = listArguments[["xfunction"]]; listArguments[["xfunction"]]=NULL #delete from the list of arguments

xsetRdataOutput = paste(thefunction,"RData",sep=".")
if (!is.null(listArguments[["xsetRdataOutput"]])){
  xsetRdataOutput = listArguments[["xsetRdataOutput"]]; listArguments[["xsetRdataOutput"]]=NULL
}

#saving the specific parameters
rplotspdf = "Rplots.pdf"
if (!is.null(listArguments[["rplotspdf"]])){
  rplotspdf = listArguments[["rplotspdf"]]; listArguments[["rplotspdf"]]=NULL
}
sampleMetadataOutput = "sampleMetadata.tsv"
if (!is.null(listArguments[["sampleMetadataOutput"]])){
  sampleMetadataOutput = listArguments[["sampleMetadataOutput"]]; listArguments[["sampleMetadataOutput"]]=NULL
}
variableMetadataOutput = "variableMetadata.tsv"
if (!is.null(listArguments[["variableMetadataOutput"]])){
  variableMetadataOutput = listArguments[["variableMetadataOutput"]]; listArguments[["variableMetadataOutput"]]=NULL
}
dataMatrixOutput = "dataMatrix.tsv"
if (!is.null(listArguments[["dataMatrixOutput"]])){
  dataMatrixOutput = listArguments[["dataMatrixOutput"]]; listArguments[["dataMatrixOutput"]]=NULL
}
if (!is.null(listArguments[["convertRTMinute"]])){
  convertRTMinute = listArguments[["convertRTMinute"]]; listArguments[["convertRTMinute"]]=NULL
}
if (!is.null(listArguments[["numDigitsMZ"]])){
  numDigitsMZ = listArguments[["numDigitsMZ"]]; listArguments[["numDigitsMZ"]]=NULL
}
if (!is.null(listArguments[["numDigitsRT"]])){
  numDigitsRT = listArguments[["numDigitsRT"]]; listArguments[["numDigitsRT"]]=NULL
}
if (!is.null(listArguments[["intval"]])){
  intval = listArguments[["intval"]]; listArguments[["intval"]]=NULL
}

if (thefunction %in% c("xcmsSet","retcor")) {
  ticspdf = listArguments[["ticspdf"]]; listArguments[["ticspdf"]]=NULL
  bicspdf = listArguments[["bicspdf"]]; listArguments[["bicspdf"]]=NULL
}

#necessary to unzip .zip file uploaded to Galaxy
#thanks to .zip file it's possible to upload many file as the same time conserving the tree hierarchy of directories


if (!is.null(listArguments[["zipfile"]])){
  zipfile= listArguments[["zipfile"]]; listArguments[["zipfile"]]=NULL
}

if (!is.null(listArguments[["singlefile_galaxyPath"]])){
  singlefile = list()
  singlefile[[listArguments[["singlefile_sampleName"]]]] = listArguments[["singlefile_galaxyPath"]]; listArguments[["singlefile_galaxyPath"]]=NULL; listArguments[["singlefile_sampleName"]]=NULL
}

if (!is.null(listArguments[["library"]])){
  directory=listArguments[["library"]]; listArguments[["library"]]=NULL
  if(!file.exists(directory)){
    error_message=paste("Cannot access the directory:",directory,". Please verify if the directory exists or not.")
    print(error_message)
    stop(error_message)
  }
}

# We unzip automatically the chromatograms from the zip files.
if (thefunction %in% c("xcmsSet","retcor","fillPeaks"))  {
  if(exists("singlefile_galaxyPath") && (singlefile_galaxyPath!="")) {
    if(!file.exists(singlefile_galaxyPath)){
      error_message=paste("Cannot access the sample:",singlefile_sampleName,"located:",singlefile_galaxyPath,". Please, contact your administrator ... if you have one!")
      print(error_message); stop(error_message)
    }

    file.symlink(singlefile_galaxyPath,singlefile_sampleName)

    directory = "."
    
    md5sumList=list("origin"=getMd5sum(directory))

  }
  if(exists("zipfile") && (zipfile!="")) {
    if(!file.exists(zipfile)){
      error_message=paste("Cannot access the Zip file:",zipfile,". Please, contact your administrator ... if you have one!")
      print(error_message)
      stop(error_message)
    }

    #list all file in the zip file
    #zip_files=unzip(zipfile,list=T)[,"Name"]

    #unzip
    suppressWarnings(unzip(zipfile, unzip="unzip"))

    #get the directory name
    filesInZip=unzip(zipfile, list=T);
    directories=unique(unlist(lapply(strsplit(filesInZip$Name,"/"), function(x) x[1])));
    directories=directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
    directory = "."
    if (length(directories) == 1) directory = directories

    cat("files_root_directory\t",directory,"\n")

    md5sumList=list("origin"=getMd5sum(directory))
  }
}

#addition of the directory to the list of arguments in the first position
if (thefunction == "xcmsSet") {
  checkXmlStructure(directory)
  checkFilesCompatibilityWithXcms(directory)
  listArguments=append(directory, listArguments)
}


#addition of xset object to the list of arguments in the first position
if (exists("xset")){
  listArguments=append(list(xset), listArguments)
}

cat("\n\n")






# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


#Verification of a group step before doing the fillpeaks job.

if (thefunction == "fillPeaks") {
  res=try(is.null(groupnames(xset)))
  if (class(res) == "try-error"){
    error<-geterrmessage()
    write(error, stderr())
    stop("You must always do a group step after a retcor. Otherwise it won't work for the fillpeaks step")
  }

}

#change the default display settings
#dev.new(file="Rplots.pdf", width=16, height=12)
pdf(file=rplotspdf, width=16, height=12)
if (thefunction == "group") {
  par(mfrow=c(2,2))
}
#else if (thefunction == "retcor") {
#try to change the legend display
#     par(xpd=NA)
#     par(xpd=T, mar=par()$mar+c(0,0,0,4))
#}


#execution of the function "thefunction" with the parameters given in "listArguments"

cat("\t\tCOMPUTE\n")
xset = do.call(thefunction, listArguments)


cat("\n\n")

dev.off() #dev.new(file="Rplots.pdf", width=16, height=12)

if (thefunction  == "xcmsSet") {

  #transform the files absolute pathways into relative pathways
  xset@filepaths<-sub(paste(getwd(),"/",sep="") ,"", xset@filepaths)

  if(exists("zipfile") && (zipfile!="")) {

    #Modify the samples names (erase the path)
    for(i in 1:length(sampnames(xset))){

      sample_name=unlist(strsplit(sampnames(xset)[i], "/"))
      sample_name=sample_name[length(sample_name)]
      sample_name= unlist(strsplit(sample_name,"[.]"))[1]
      sampnames(xset)[i]=sample_name

    }

  }

}

# -- TIC --
if (thefunction == "xcmsSet") {
  cat("\t\tGET TIC GRAPH\n")
  sampleNamesList = getSampleMetadata(xcmsSet=xset, sampleMetadataOutput=sampleMetadataOutput)
  getTICs(xcmsSet=xset, pdfname=ticspdf,rt="raw")
  getBPCs(xcmsSet=xset,rt="raw",pdfname=bicspdf)
} else if (thefunction == "retcor") {
  cat("\t\tGET TIC GRAPH\n")
  getTICs(xcmsSet=xset, pdfname=ticspdf,rt="corrected")
  getBPCs(xcmsSet=xset,rt="corrected",pdfname=bicspdf)
}

if (thefunction == "fillPeaks" && exists("intval")) {
  cat("\t\tGET THE PEAK LIST\n")
  getPeaklistW4M(xset,intval,convertRTMinute,numDigitsMZ,numDigitsRT,variableMetadataOutput,dataMatrixOutput)
}


cat("\n\n")

# ----- EXPORT -----

cat("\tXSET OBJECT INFO\n")
print(xset)
#delete the parameters to avoid the passage to the next tool in .RData image


#saving R data in .Rdata file to save the variables used in the present tool
objects2save = c("xset","zipfile","singlefile_galaxyPath","singlefile_sampleName","listOFlistArguments","md5sumList","sampleNamesList")
save(list=objects2save[objects2save %in% ls()], file=xsetRdataOutput)

cat("\n\n")


cat("\tDONE\n")

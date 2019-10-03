#!/usr/bin/env Rscript
#Authors Gildas Le Corguille and Yann Guitton


# ----- LOG FILE -----
log_file=file("log.txt", open = "wt")
sink(log_file)
sink(log_file, type = "output")


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname){ argv <- commandArgs(trailingOnly=FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep="/")) }
source_local("lib.r")

pkgs <- c("IPO","batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");


# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n")
args = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(args), col.names=F, quote=F, sep='\t')

cat("\n\n");


# ----- ARGUMENTS PROCESSING -----
cat("\tINFILE PROCESSING INFO\n")
options(bitmapType='cairo')

#image is an .RData file necessary to use xset variable given by previous tools
if (!is.null(args$image)){
  load(args$image); args$image=NULL
}

cat("\n\n")

#Import the different functions

# ----- PROCESSING INFILE -----
cat("\tARGUMENTS PROCESSING INFO\n")


parametersOutput = "parametersOutput.tsv"
if (!is.null(args$parametersOutput)){
  parametersOutput = args$parametersOutput; args$parametersOutput=NULL
}

samplebyclass = 2
if (!is.null(args$samplebyclass)){
  samplebyclass = args$samplebyclass; args$samplebyclass=NULL
}

#necessary to unzip .zip file uploaded to Galaxy
#thanks to .zip file it's possible to upload many file as the same time conserving the tree hierarchy of directories


if (!is.null(args$zipfile)){
  zipfile= args$zipfile; args$zipfile=NULL
}


if (!is.null(args$singlefile_galaxyPath)){
    singlefile_galaxyPath = unlist(strsplit(args$singlefile_galaxyPath,",")); args$singlefile_galaxyPath=NULL
    singlefile_sampleName = unlist(strsplit(args$singlefile_sampleName,",")); args$singlefile_sampleName=NULL
}

# single file case
#@TODO: need to be refactoring
if(exists("singlefile_galaxyPath") && (singlefile_galaxyPath!="")) {

    cwd=getwd()
    dir.create("raw")
    setwd("raw")

    for (singlefile_galaxyPath_i in seq(1:length(singlefile_galaxyPath))) {
        if(!file.exists(singlefile_galaxyPath[singlefile_galaxyPath_i])){
            error_message=paste("Cannot access the sample:",singlefile_sampleName[singlefile_galaxyPath_i],"located:",singlefile_galaxyPath[singlefile_galaxyPath_i],". Please, contact your administrator ... if you have one!")
            print(error_message); stop(error_message)
        }
        file.symlink(singlefile_galaxyPath[singlefile_galaxyPath_i],singlefile_sampleName[singlefile_galaxyPath_i])
    }

    setwd(cwd)

    directory = "raw"

}

# We unzip automatically the chromatograms from the zip files.
if(exists("zipfile") && (zipfile!="")) {
    if(!file.exists(zipfile)){
        error_message=paste("Cannot access the Zip file:",zipfile,". Please, contact your administrator ... if you have one!")
        print(error_message)
        stop(error_message)
    }

    #list all file in the zip file
    #zip_files=unzip(zipfile,list=T)[,"Name"]

    # Because IPO only want raw data in its working directory
    dir.create("ipoworkingdir")
    setwd("ipoworkingdir")

    #unzip
    suppressWarnings(unzip(zipfile, unzip="unzip"))

    #get the directory name
    filesInZip=unzip(zipfile, list=T);
    directories=unique(unlist(lapply(strsplit(filesInZip$Name,"/"), function(x) x[1])));
    directories=directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
    directory = "."
    if (length(directories) == 1) directory = directories

    cat("files_root_directory\t",directory,"\n")


}

#addition of the directory to the list of arguments in the first position
checkXmlStructure(directory)
checkFilesCompatibilityWithXcms(directory)

cat("\n\n")






# ----- MAIN PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")


ipo4retgroup(xset, directory, parametersOutput, args, samplebyclass)



cat("\n\n")


# ----- EXPORT -----

#cat("\tEXPORTING INFO\n")

#save.image(file="ipo-retcor.RData")

#cat("\n\n")


cat("\tDONE\n")

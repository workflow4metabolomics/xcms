#@authors ABiMS TEAM, Y. Guitton
# lib.r for Galaxy Workflow4Metabolomics xcms tools

#@author G. Le Corguille
# solve an issue with batch if arguments are logical TRUE/FALSE
parseCommandArgs <- function(...) {
    args <- batch::parseCommandArgs(...)
    for (key in names(args)) {
        if (args[key] %in% c("TRUE","FALSE"))
            args[key] = as.logical(args[key])
    }
    return(args)
}

#@author G. Le Corguille
# This function will
# - load the packages
# - display the sessionInfo
loadAndDisplayPackages <- function(pkgs) {
    for(pkg in pkgs) suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))

    sessioninfo = sessionInfo()
    cat(sessioninfo$R.version$version.string,"\n")
    cat("Main packages:\n")
    for (pkg in names(sessioninfo$otherPkgs)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
    cat("Other loaded packages:\n")
    for (pkg in names(sessioninfo$loadedOnly)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
}

#@author G. Le Corguille
# This function merge several chromBPI or chromTIC into one.
mergeChrom <- function(chrom_merged, chrom) {
    if (is.null(chrom_merged)) return(NULL)
    chrom_merged@.Data <- cbind(chrom_merged@.Data, chrom@.Data)
    return(chrom_merged)
}

#@author G. Le Corguille
# This function merge several xdata into one.
mergeXData <- function(args) {
    chromTIC <- NULL
    chromBPI <- NULL
    chromTIC_adjusted <- NULL
    chromBPI_adjusted <- NULL
    for(image in args$images) {

        load(image)
        # Handle infiles
        if (!exists("singlefile")) singlefile <- NULL
        if (!exists("zipfile")) zipfile <- NULL
        rawFilePath <- getRawfilePathFromArguments(singlefile, zipfile, args)
        zipfile <- rawFilePath$zipfile
        singlefile <- rawFilePath$singlefile
        retrieveRawfileInTheWorkingDirectory(singlefile, zipfile)

        if (exists("raw_data")) xdata <- raw_data
        if (!exists("xdata")) stop("\n\nERROR: The RData doesn't contain any object called 'xdata'. This RData should have been created by an old version of XMCS 2.*")

        cat(sampleNamesList$sampleNamesOrigin,"\n")

        if (!exists("xdata_merged")) {
            xdata_merged <- xdata
            singlefile_merged <- singlefile
            md5sumList_merged <- md5sumList
            sampleNamesList_merged <- sampleNamesList
            chromTIC_merged <- chromTIC
            chromBPI_merged <- chromBPI
            chromTIC_adjusted_merged <- chromTIC_adjusted
            chromBPI_adjusted_merged <- chromBPI_adjusted
        } else {
            if (is(xdata, "XCMSnExp")) xdata_merged <- c(xdata_merged,xdata)
            else if (is(xdata, "OnDiskMSnExp")) xdata_merged <- xcms:::.concatenate_OnDiskMSnExp(xdata_merged,xdata)
            else stop("\n\nERROR: The RData either a OnDiskMSnExp object called raw_data or a XCMSnExp object called xdata")

            singlefile_merged <- c(singlefile_merged,singlefile)
            md5sumList_merged$origin <- rbind(md5sumList_merged$origin,md5sumList$origin)
            sampleNamesList_merged$sampleNamesOrigin <- c(sampleNamesList_merged$sampleNamesOrigin,sampleNamesList$sampleNamesOrigin)
            sampleNamesList_merged$sampleNamesMakeNames <- c(sampleNamesList_merged$sampleNamesMakeNames,sampleNamesList$sampleNamesMakeNames)
            chromTIC_merged <- mergeChrom(chromTIC_merged, chromTIC)
            chromBPI_merged <- mergeChrom(chromBPI_merged, chromBPI)
            chromTIC_adjusted_merged <- mergeChrom(chromTIC_adjusted_merged, chromTIC_adjusted)
            chromBPI_adjusted_merged <- mergeChrom(chromBPI_adjusted_merged, chromBPI_adjusted)
        }
    }
    rm(image)
    xdata <- xdata_merged; rm(xdata_merged)
    singlefile <- singlefile_merged; rm(singlefile_merged)
    md5sumList <- md5sumList_merged; rm(md5sumList_merged)
    sampleNamesList <- sampleNamesList_merged; rm(sampleNamesList_merged)

    if (!is.null(args$sampleMetadata)) {
        cat("\tXSET PHENODATA SETTING...\n")
        sampleMetadataFile <- args$sampleMetadata
        sampleMetadata <- getDataFrameFromFile(sampleMetadataFile, header=F)
        xdata@phenoData@data$sample_group=sampleMetadata$V2[match(xdata@phenoData@data$sample_name,sampleMetadata$V1)]

        if (any(is.na(pData(xdata)$sample_group))) {
            sample_missing <- pData(xdata)$sample_name[is.na(pData(xdata)$sample_group)]
            error_message <- paste("Those samples are missing in your sampleMetadata:", paste(sample_missing, collapse=" "))
            print(error_message)
            stop(error_message)
        }
    }

    if (!is.null(chromTIC_merged)) { chromTIC <- chromTIC_merged; chromTIC@phenoData <- xdata@phenoData }
    if (!is.null(chromBPI_merged)) { chromBPI <- chromBPI_merged; chromBPI@phenoData <- xdata@phenoData }
    if (!is.null(chromTIC_adjusted_merged)) { chromTIC_adjusted <- chromTIC_adjusted_merged; chromTIC_adjusted@phenoData <- xdata@phenoData }
    if (!is.null(chromBPI_adjusted_merged)) { chromBPI_adjusted <- chromBPI_adjusted_merged; chromBPI_adjusted@phenoData <- xdata@phenoData }

    return(list("xdata"=xdata, "singlefile"=singlefile, "md5sumList"=md5sumList,"sampleNamesList"=sampleNamesList, "chromTIC"=chromTIC, "chromBPI"=chromBPI, "chromTIC_adjusted"=chromTIC_adjusted, "chromBPI_adjusted"=chromBPI_adjusted))
}

#@author G. Le Corguille
# This function convert if it is required the Retention Time in minutes
RTSecondToMinute <- function(variableMetadata, convertRTMinute) {
    if (convertRTMinute){
        #converting the retention times (seconds) into minutes
        print("converting the retention times into minutes in the variableMetadata")
        variableMetadata[,"rt"] <- variableMetadata[,"rt"]/60
        variableMetadata[,"rtmin"] <- variableMetadata[,"rtmin"]/60
        variableMetadata[,"rtmax"] <- variableMetadata[,"rtmax"]/60
    }
    return (variableMetadata)
}

#@author G. Le Corguille
# This function format ions identifiers
formatIonIdentifiers <- function(variableMetadata, numDigitsRT=0, numDigitsMZ=0) {
    splitDeco <- strsplit(as.character(variableMetadata$name),"_")
    idsDeco <- sapply(splitDeco, function(x) { deco=unlist(x)[2]; if (is.na(deco)) return ("") else return(paste0("_",deco)) })
    namecustom <- make.unique(paste0("M",round(variableMetadata[,"mz"],numDigitsMZ),"T",round(variableMetadata[,"rt"],numDigitsRT),idsDeco))
    variableMetadata <- cbind(name=variableMetadata$name, namecustom=namecustom, variableMetadata[,!(colnames(variableMetadata) %in% c("name"))])
    return(variableMetadata)
}

#@author G. Le Corguille
# This function convert the remain NA to 0 in the dataMatrix
naTOzeroDataMatrix <- function(dataMatrix, naTOzero) {
    if (naTOzero){
        dataMatrix[is.na(dataMatrix)] <- 0
    }
    return (dataMatrix)
}

#@author G. Le Corguille
# Draw the plotChromPeakDensity 3 per page in a pdf file
getPlotChromPeakDensity <- function(xdata, param = NULL, mzdigit=4) {
    pdf(file="plotChromPeakDensity.pdf", width=16, height=12)

    par(mfrow = c(3, 1), mar = c(4, 4, 1, 0.5))

    group_colors <- brewer.pal(3, "Set1")[1:length(unique(xdata$sample_group))]
    names(group_colors) <- unique(xdata$sample_group)

    xlim <- c(min(featureDefinitions(xdata)$rtmin), max(featureDefinitions(xdata)$rtmax))
    for (i in 1:nrow(featureDefinitions(xdata))) {
        mzmin = featureDefinitions(xdata)[i,]$mzmin
        mzmax = featureDefinitions(xdata)[i,]$mzmax
        plotChromPeakDensity(xdata, param = param, mz=c(mzmin,mzmax), col=group_colors, pch=16, xlim=xlim, main=paste(round(mzmin,mzdigit),round(mzmax,mzdigit)))
        legend("topright", legend=names(group_colors), col=group_colors, cex=0.8, lty=1)
    }

    dev.off()
}

#@author G. Le Corguille
# Draw the plotChromPeakDensity 3 per page in a pdf file
getPlotAdjustedRtime <- function(xdata) {

    pdf(file="raw_vs_adjusted_rt.pdf", width=16, height=12)

    # Color by group
    group_colors <- brewer.pal(3, "Set1")[1:length(unique(xdata$sample_group))]
    if (length(group_colors) > 1) {
        names(group_colors) <- unique(xdata$sample_group)
        plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group])
        legend("topright", legend=names(group_colors), col=group_colors, cex=0.8, lty=1)
    }

    # Color by sample
    plotAdjustedRtime(xdata, col = rainbow(length(xdata@phenoData@data$sample_name)))
    legend("topright", legend=xdata@phenoData@data$sample_name, col=rainbow(length(xdata@phenoData@data$sample_name)), cex=0.8, lty=1)

    dev.off()
}

#@author G. Le Corguille
# value: intensity values to be used into, maxo or intb
getPeaklistW4M <- function(xdata, intval="into", convertRTMinute=F, numDigitsMZ=4, numDigitsRT=0, naTOzero=T, variableMetadataOutput, dataMatrixOutput, sampleNamesList) {
    dataMatrix <- featureValues(xdata, method="medret", value=intval)
    colnames(dataMatrix) <- make.names(tools::file_path_sans_ext(colnames(dataMatrix)))
    dataMatrix = cbind(name=groupnames(xdata), dataMatrix)
    variableMetadata <- featureDefinitions(xdata)
    colnames(variableMetadata)[1] = "mz"; colnames(variableMetadata)[4] = "rt"
    variableMetadata = data.frame(name=groupnames(xdata), variableMetadata)

    variableMetadata <- RTSecondToMinute(variableMetadata, convertRTMinute)
    variableMetadata <- formatIonIdentifiers(variableMetadata, numDigitsRT=numDigitsRT, numDigitsMZ=numDigitsMZ)
    dataMatrix <- naTOzeroDataMatrix(dataMatrix, naTOzero)

    # FIX: issue when the vector at peakidx is too long and is written in a new line during the export
    variableMetadata[,"peakidx"] <- vapply(variableMetadata[,"peakidx"], FUN = paste, FUN.VALUE = character(1), collapse = ",")

    write.table(variableMetadata, file=variableMetadataOutput,sep="\t",quote=F,row.names=F)
    write.table(dataMatrix, file=dataMatrixOutput,sep="\t",quote=F,row.names=F)

}

#@author G. Le Corguille
# It allow different of field separators
getDataFrameFromFile <- function(filename, header=T) {
    myDataFrame <- read.table(filename, header=header, sep=";", stringsAsFactors=F)
    if (ncol(myDataFrame) < 2) myDataFrame <- read.table(filename, header=header, sep="\t", stringsAsFactors=F)
    if (ncol(myDataFrame) < 2) myDataFrame <- read.table(filename, header=header, sep=",", stringsAsFactors=F)
    if (ncol(myDataFrame) < 2) {
        error_message="Your tabular file seems not well formatted. The column separators accepted are ; , and tabulation"
        print(error_message)
        stop(error_message)
    }
    return(myDataFrame)
}

#@author G. Le Corguille
# Draw the BPI and TIC graphics
# colored by sample names or class names
getPlotChromatogram <- function(chrom, xdata, pdfname="Chromatogram.pdf", aggregationFun = "max") {

    if (aggregationFun == "sum")
        type="Total Ion Chromatograms"
    else
        type="Base Peak Intensity Chromatograms"

    adjusted="Raw"
    if (hasAdjustedRtime(xdata))
        adjusted="Adjusted"

    main <- paste(type,":",adjusted,"data")

    pdf(pdfname, width=16, height=10)

    # Color by group
    group_colors <- brewer.pal(3, "Set1")[1:length(unique(xdata$sample_group))]
    if (length(group_colors) > 1) {
        names(group_colors) <- unique(xdata$sample_group)
        plot(chrom, col = group_colors[chrom$sample_group], main=main)
        legend("topright", legend=names(group_colors), col=group_colors, cex=0.8, lty=1)
    }

    # Color by sample
    plot(chrom, col = rainbow(length(xdata@phenoData@data$sample_name)), main=main)
    legend("topright", legend=xdata@phenoData@data$sample_name, col=rainbow(length(xdata@phenoData@data$sample_name)), cex=0.8, lty=1)

    dev.off()
}


# Get the polarities from all the samples of a condition
#@author Misharl Monsoor misharl.monsoor@sb-roscoff.fr ABiMS TEAM
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr ABiMS TEAM
getSampleMetadata <- function(xdata=NULL, sampleMetadataOutput="sampleMetadata.tsv") {
    cat("Creating the sampleMetadata file...\n")

    #Create the sampleMetada dataframe
    sampleMetadata <- xdata@phenoData@data
    rownames(sampleMetadata) <- NULL
    colnames(sampleMetadata) <-  c("sample_name", "class")

    sampleNamesOrigin <- sampleMetadata$sample_name
    sampleNamesMakeNames <- make.names(sampleNamesOrigin)

    if (any(duplicated(sampleNamesMakeNames))) {
        write("\n\nERROR: Usually, R has trouble to deal with special characters in its column names, so it rename them using make.names().\nIn your case, at least two columns after the renaming obtain the same name, thus XCMS will collapse those columns per name.", stderr())
        for (sampleName in sampleNamesOrigin) {
            write(paste(sampleName,"\t->\t",make.names(sampleName)),stderr())
        }
        stop("\n\nERROR: One or more of your files will not be import by xcmsSet. It may due to bad characters in their filenames.")
    }

    if (!all(sampleNamesOrigin == sampleNamesMakeNames)) {
        cat("\n\nWARNING: Usually, R has trouble to deal with special characters in its column names, so it rename them using make.names()\nIn your case, one or more sample names will be renamed in the sampleMetadata and dataMatrix files:\n")
        for (sampleName in sampleNamesOrigin) {
            cat(paste(sampleName,"\t->\t",make.names(sampleName),"\n"))
        }
    }

    sampleMetadata$sample_name <- sampleNamesMakeNames


    #For each sample file, the following actions are done
    for (fileIdx in 1:length(fileNames(xdata))) {
        #Check if the file is in the CDF format
        if (!mzR:::netCDFIsFile(fileNames(xdata))) {

            # If the column isn't exist, with add one filled with NA
            if (is.null(sampleMetadata$polarity)) sampleMetadata$polarity <- NA

            #Extract the polarity (a list of polarities)
            polarity <- fData(xdata)[fData(xdata)$fileIdx == fileIdx,"polarity"]
            #Verify if all the scans have the same polarity
            uniq_list <- unique(polarity)
            if (length(uniq_list)>1){
                polarity <- "mixed"
            } else {
                polarity <- as.character(uniq_list)
            }

            #Set the polarity attribute
            sampleMetadata$polarity[fileIdx] <- polarity
        }

    }

    write.table(sampleMetadata, sep="\t", quote=FALSE, row.names=FALSE, file=sampleMetadataOutput)

    return(list("sampleNamesOrigin"=sampleNamesOrigin, "sampleNamesMakeNames"=sampleNamesMakeNames))

}


# This function check if xcms will found all the files
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr ABiMS TEAM
checkFilesCompatibilityWithXcms <- function(directory) {
    cat("Checking files filenames compatibilities with xmcs...\n")
    # WHAT XCMS WILL FIND
    filepattern <- c("([Nn][Ee][Tt])?[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]","[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep=""),collapse="|")
    info <- file.info(directory)
    listed <- list.files(directory[info$isdir], pattern=filepattern, recursive=TRUE, full.names=TRUE)
    files <- c(directory[!info$isdir], listed)
    files_abs <- file.path(getwd(), files)
    exists <- file.exists(files_abs)
    files[exists] <- files_abs[exists]
    files[exists] <- sub("//","/",files[exists])

    # WHAT IS ON THE FILESYSTEM
    filesystem_filepaths <- system(paste0("find \"",getwd(),"/",directory,"\" -not -name '\\.*' -not -path '*conda-env*' -type f -name \"*\""), intern=T)
    filesystem_filepaths <- filesystem_filepaths[grep(filepattern, filesystem_filepaths, perl=T)]

    # COMPARISON
    if (!is.na(table(filesystem_filepaths %in% files)["FALSE"])) {
        write("\n\nERROR: List of the files which will not be imported by xcmsSet",stderr())
        write(filesystem_filepaths[!(filesystem_filepaths %in% files)],stderr())
        stop("\n\nERROR: One or more of your files will not be import by xcmsSet. It may due to bad characters in their filenames.")
    }
}


#This function list the compatible files within the directory as xcms did
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr ABiMS TEAM
getMSFiles <- function (directory) {
    filepattern <- c("([Nn][Ee][Tt])?[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]","[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep=""),collapse="|")
    info <- file.info(directory)
    listed <- list.files(directory[info$isdir], pattern=filepattern,recursive=TRUE, full.names=TRUE)
    files <- c(directory[!info$isdir], listed)
    exists <- file.exists(files)
    files <- files[exists]
    return(files)
}

# This function check if XML contains special caracters. It also checks integrity and completness.
#@author Misharl Monsoor misharl.monsoor@sb-roscoff.fr ABiMS TEAM
checkXmlStructure <- function (directory) {
    cat("Checking XML structure...\n")

    cmd <- paste0("IFS=$'\n'; for xml in $(find '",directory,"' -not -name '\\.*' -not -path '*conda-env*' -type f -iname '*.*ml*'); do if [ $(xmllint --nonet --noout \"$xml\" 2> /dev/null; echo $?) -gt 0 ]; then echo $xml;fi; done;")
    capture <- system(cmd, intern=TRUE)

    if (length(capture)>0){
        #message=paste("The following mzXML or mzML file is incorrect, please check these files first:",capture)
        write("\n\nERROR: The following mzXML or mzML file(s) are incorrect, please check these files first:", stderr())
        write(capture, stderr())
        stop("ERROR: xcmsSet cannot continue with incorrect mzXML or mzML files")
    }

}


# This function check if XML contain special characters
#@author Misharl Monsoor misharl.monsoor@sb-roscoff.fr ABiMS TEAM
deleteXmlBadCharacters<- function (directory) {
    cat("Checking Non ASCII characters in the XML...\n")

    processed <- F
    l <- system( paste0("find '",directory, "' -not -name '\\.*' -not -path '*conda-env*' -type f -iname '*.*ml*'"), intern=TRUE)
    for (i in l){
        cmd <- paste("LC_ALL=C grep '[^ -~]' \"", i, "\"", sep="")
        capture <- suppressWarnings(system(cmd, intern=TRUE))
        if (length(capture)>0){
            cmd <- paste("perl -i -pe 's/[^[:ascii:]]//g;'",i)
            print( paste("WARNING: Non ASCII characters have been removed from the ",i,"file") )
            c <- system(cmd, intern=TRUE)
            capture <- ""
            processed <- T
        }
    }
    if (processed) cat("\n\n")
    return(processed)
}


# This function will compute MD5 checksum to check the data integrity
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getMd5sum <- function (directory) {
    cat("Compute md5 checksum...\n")
    # WHAT XCMS WILL FIND
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]","[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep=""),collapse="|")
    info <- file.info(directory)
    listed <- list.files(directory[info$isdir], pattern=filepattern, recursive=TRUE, full.names=TRUE)
    files <- c(directory[!info$isdir], listed)
    exists <- file.exists(files)
    files <- files[exists]

    library(tools)

    #cat("\n\n")

    return(as.matrix(md5sum(files)))
}


# This function get the raw file path from the arguments
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getRawfilePathFromArguments <- function(singlefile, zipfile, args, prefix="") {
  if (!(prefix %in% c("","Positive","Negative","MS1","MS2"))) stop("prefix must be either '', 'Positive', 'Negative', 'MS1' or 'MS2'")

  if (!is.null(args[[paste0("zipfile",prefix)]])) zipfile <- args[[paste0("zipfile",prefix)]]

  if (!is.null(args[[paste0("singlefile_galaxyPath",prefix)]])) {
    singlefile_galaxyPaths <- args[[paste0("singlefile_galaxyPath",prefix)]]
    singlefile_sampleNames <- args[[paste0("singlefile_sampleName",prefix)]]
  }
  if (exists("singlefile_galaxyPaths")){
    singlefile_galaxyPaths <- unlist(strsplit(singlefile_galaxyPaths,"\\|"))
    singlefile_sampleNames <- unlist(strsplit(singlefile_sampleNames,"\\|"))

    singlefile <- NULL
    for (singlefile_galaxyPath_i in seq(1:length(singlefile_galaxyPaths))) {
      singlefile_galaxyPath <- singlefile_galaxyPaths[singlefile_galaxyPath_i]
      singlefile_sampleName <- singlefile_sampleNames[singlefile_galaxyPath_i]
      # In case, an url is used to import data within Galaxy
      singlefile_sampleName <- tail(unlist(strsplit(singlefile_sampleName,"/")), n=1)
      singlefile[[singlefile_sampleName]] <- singlefile_galaxyPath
    }
  }
  return(list(zipfile=zipfile, singlefile=singlefile))
}

# This function retrieve the raw file in the working directory
#   - if zipfile: unzip the file with its directory tree
#   - if singlefiles: set symlink with the good filename
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
retrieveRawfileInTheWorkingDirectory <- function(singlefile, zipfile) {
    if(!is.null(singlefile) && (length("singlefile")>0)) {
        for (singlefile_sampleName in names(singlefile)) {
            singlefile_galaxyPath <- singlefile[[singlefile_sampleName]]
            if(!file.exists(singlefile_galaxyPath)){
                error_message <- paste("Cannot access the sample:",singlefile_sampleName,"located:",singlefile_galaxyPath,". Please, contact your administrator ... if you have one!")
                print(error_message); stop(error_message)
            }

            if (!suppressWarnings( try (file.link(singlefile_galaxyPath, singlefile_sampleName), silent=T)))
                file.copy(singlefile_galaxyPath, singlefile_sampleName)

        }
        directory <- "."

    }
    if(!is.null(zipfile) && (zipfile != "")) {
        if(!file.exists(zipfile)){
            error_message <- paste("Cannot access the Zip file:",zipfile,". Please, contact your administrator ... if you have one!")
            print(error_message)
            stop(error_message)
        }

        #list all file in the zip file
        #zip_files <- unzip(zipfile,list=T)[,"Name"]

        #unzip
        suppressWarnings(unzip(zipfile, unzip="unzip"))

        #get the directory name
        suppressWarnings(filesInZip <- unzip(zipfile, list=T))
        directories <- unique(unlist(lapply(strsplit(filesInZip$Name,"/"), function(x) x[1])))
        directories <- directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
        directory <- "."
        if (length(directories) == 1) directory <- directories

        cat("files_root_directory\t",directory,"\n")

    }
    return (directory)
}


# This function retrieve a xset like object
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getxcmsSetObject <- function(xobject) {
    # XCMS 1.x
    if (class(xobject) == "xcmsSet")
        return (xobject)
    # XCMS 3.x
    if (class(xobject) == "XCMSnExp") {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xobject, 'xcmsSet'))
        if (!is.null(xset@phenoData$sample_group))
            sampclass(xset) <- xset@phenoData$sample_group
        else
            sampclass(xset) <- "."
        return (xset)
    }
}

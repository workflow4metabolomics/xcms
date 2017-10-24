#Authors ABiMS TEAM
#Lib.r for Galaxy Workflow4Metabolomics xcms tools
#
#version 2.4: lecorguille
#   add getPeaklistW4M
#version 2.3: yguitton
#   correction for empty PDF when only 1 class
#version 2.2
#   correct bug in Base Peak Chromatogram (BPC) option, not only TIC when scanrange used in xcmsSet
#   Note if scanrange is used a warning is prompted in R console but do not stop PDF generation
#version 2.1: yguitton
#   Modifications made by Guitton Yann


#@author G. Le Corguille
#This function convert if it is required the Retention Time in minutes
RTSecondToMinute <- function(variableMetadata, convertRTMinute) {
    if (convertRTMinute){
        #converting the retention times (seconds) into minutes
        print("converting the retention times into minutes in the variableMetadata")
        variableMetadata[,"rt"]=variableMetadata[,"rt"]/60
        variableMetadata[,"rtmin"]=variableMetadata[,"rtmin"]/60
        variableMetadata[,"rtmax"]=variableMetadata[,"rtmax"]/60
    }
    return (variableMetadata)
}

#@author G. Le Corguille
#This function format ions identifiers
formatIonIdentifiers <- function(variableMetadata, numDigitsRT=0, numDigitsMZ=0) {
    splitDeco = strsplit(as.character(variableMetadata$name),"_")
    idsDeco = sapply(splitDeco, function(x) { deco=unlist(x)[2]; if (is.na(deco)) return ("") else return(paste0("_",deco)) })
    namecustom = make.unique(paste0("M",round(variableMetadata[,"mz"],numDigitsMZ),"T",round(variableMetadata[,"rt"],numDigitsRT),idsDeco))
    variableMetadata=cbind(name=variableMetadata$name, namecustom=namecustom, variableMetadata[,!(colnames(variableMetadata) %in% c("name"))])
    return(variableMetadata)
}

#@author G. Le Corguille
# value: intensity values to be used into, maxo or intb
getPeaklistW4M <- function(xset, intval="into",convertRTMinute=F,numDigitsMZ=4,numDigitsRT=0,variableMetadataOutput,dataMatrixOutput) {
    variableMetadata_dataMatrix = peakTable(xset, method="medret", value=intval)
    variableMetadata_dataMatrix = cbind(name=groupnames(xset),variableMetadata_dataMatrix)

    dataMatrix = variableMetadata_dataMatrix[,(make.names(colnames(variableMetadata_dataMatrix)) %in% c("name", make.names(sampnames(xset))))]

    variableMetadata = variableMetadata_dataMatrix[,!(make.names(colnames(variableMetadata_dataMatrix)) %in% c(make.names(sampnames(xset))))]
    variableMetadata = RTSecondToMinute(variableMetadata, convertRTMinute)
    variableMetadata = formatIonIdentifiers(variableMetadata, numDigitsRT=numDigitsRT, numDigitsMZ=numDigitsMZ)

    write.table(variableMetadata, file=variableMetadataOutput,sep="\t",quote=F,row.names=F)
    write.table(dataMatrix, file=dataMatrixOutput,sep="\t",quote=F,row.names=F)
}

#@author G. Le Corguille
exportTicBpcTabular <- function(dataset, filenameBase, ticORbpc, rt='raw') {

        rawORcorrected = ''
        title = ''
        if (rt=='corrected') {
            rawORcorrected = '_corrected'
            title = ' corrected by retcor'
        }

        if (ticORbpc == "TIC") {
            section_name = paste0('TIC',rawORcorrected)
            title = paste0('Total Ion Current (TIC) chromatogram',title)
            description = 'Sum of intensity (Y) of all ions detected at each retention time(X)'
        } else if (ticORbpc == "BPC") {
            section_name = paste0('BPC',rawORcorrected)
            title = paste0('Base Peak Chromatogram (BPC)',title)
            description = 'Sum of intensity (Y) of the most intense peaks at each retention time(X)'
        }

        filename=paste0(basename(file_path_sans_ext(filenameBase)),"-",ticORbpc,rawORcorrected,"_mqc.out")

        # Headers for MultiQC
        cat("# file_format: 'tsv'\n", sep="", file = filename)
        cat("# section_name: '",section_name,"'\n", sep="", file = filename, append = T)
        cat("# title: '",title,"'\n", sep="", file = filename, append = T)
        cat("# description: '",description,"'\n", sep="", file = filename, append = T)
        cat("# plot_type: 'linegraph'\n", sep="", file = filename, append = T)
        cat("# pconfig:\n", sep="", file = filename, append = T)
        cat("#     id: '",ticORbpc,rawORcorrected,"_lineplot'\n", sep="", file = filename, append = T)
        cat("#     ylab: 'Base Peak Intensity'\n", sep="", file = filename, append = T)
        cat("#     xlab: 'Retention Time'\n", sep="", file = filename, append = T)
        colnames(dataset) = c("Intensity","RT")
        write.table(dataset, filename ,row.names = F, sep = "\t", append = T, quote = F)
}

#@author Y. Guitton
getBPC <- function(file,rtcor=NULL, ...) {
    object <- xcmsRaw(file)
    sel <- profRange(object, ...)
    cbind(if (is.null(rtcor)) object@scantime[sel$scanidx] else rtcor ,xcms:::colMax(object@env$profile[sel$massidx,sel$scanidx,drop=FALSE]))
    #plotChrom(xcmsRaw(file), base=T)
}

#@author Y. Guitton
getBPCs <- function (xcmsSet=NULL, pdfname="BPCs.pdf",rt=c("raw","corrected"), scanrange=NULL) {
    cat("Creating BIC pdf...\n")

    if (is.null(xcmsSet)) {
        cat("Enter an xcmsSet \n")
        stop()
    } else {
        files <- filepaths(xcmsSet)
    }

    phenoDataClass<-as.vector(levels(xcmsSet@phenoData[,1])) #sometime phenoData have more than 1 column use first as class

    classnames<-vector("list",length(phenoDataClass))
    for (i in 1:length(phenoDataClass)){
        classnames[[i]]<-which( xcmsSet@phenoData[,1]==phenoDataClass[i])
    }

    N <- dim(phenoData(xcmsSet))[1]

    TIC <- vector("list",N)


    for (j in 1:N) {

        TIC[[j]] <- getBPC(files[j])
        #good for raw
        # seems strange for corrected
        #errors if scanrange used in xcmsSetgeneration
        if (!is.null(xcmsSet) && rt == "corrected")
            rtcor <- xcmsSet@rt$corrected[[j]]
        else
            rtcor <- NULL

        TIC[[j]] <- getBPC(files[j],rtcor=rtcor)

        exportTicBpcTabular(TIC[[j]], files[j], "BPC", rt=rt)

        # TIC[[j]][,1]<-rtcor
    }



    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in BPCs
    xlim = range(sapply(TIC, function(x) range(x[,1])))
    ylim = range(sapply(TIC, function(x) range(x[,2])))
    ylim = c(-ylim[2], ylim[2])


    ##plot start

    if (length(phenoDataClass)>2){
        for (k in 1:(length(phenoDataClass)-1)){
            for (l in (k+1):length(phenoDataClass)){
                #print(paste(phenoDataClass[k],"vs",phenoDataClass[l],sep=" "))
                plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",phenoDataClass[k]," vs ",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")
                colvect<-NULL
                for (j in 1:length(classnames[[k]])) {
                    tic <- TIC[[classnames[[k]][j]]]
                    # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
                    points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[k]][j]])
                }
                for (j in 1:length(classnames[[l]])) {
                    # i=class2names[j]
                    tic <- TIC[[classnames[[l]][j]]]
                    points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[l]][j]])
                }
                legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)
            }
        }
    }#end if length >2

    if (length(phenoDataClass)==2){
        k=1
        l=2
        colvect<-NULL
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",phenoDataClass[k],"vs",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")

        for (j in 1:length(classnames[[k]])) {

            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        for (j in 1:length(classnames[[l]])) {
            # i=class2names[j]
            tic <- TIC[[classnames[[l]][j]]]
            points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[l]][j]])
        }
        legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==2

    #case where only one class
    if (length(phenoDataClass)==1){
        k=1
        ylim = range(sapply(TIC, function(x) range(x[,2])))
        colvect<-NULL
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",phenoDataClass[k], sep=""), xlab = "Retention Time (min)", ylab = "BPC")

        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }

        legend("topright",paste(basename(files[c(classnames[[k]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==1

    dev.off() #pdf(pdfname,w=16,h=10)

    invisible(TIC)
}



#@author Y. Guitton
getTIC <- function(file,rtcor=NULL) {
    object <- xcmsRaw(file)
    cbind(if (is.null(rtcor)) object@scantime else rtcor, rawEIC(object,mzrange=range(object@env$mz))$intensity)
}

##
##  overlay TIC from all files in current folder or from xcmsSet, create pdf
##
#@author Y. Guitton
getTICs <- function(xcmsSet=NULL,files=NULL, pdfname="TICs.pdf",rt=c("raw","corrected")) {
    cat("Creating TIC pdf...\n")

    if (is.null(xcmsSet)) {
        filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
        filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
        if (is.null(files))
            files <- getwd()
        info <- file.info(files)
        listed <- list.files(files[info$isdir], pattern = filepattern, recursive = TRUE, full.names = TRUE)
        files <- c(files[!info$isdir], listed)
    } else {
        files <- filepaths(xcmsSet)
    }

    phenoDataClass<-as.vector(levels(xcmsSet@phenoData[,1])) #sometime phenoData have more than 1 column use first as class
    classnames<-vector("list",length(phenoDataClass))
    for (i in 1:length(phenoDataClass)){
        classnames[[i]]<-which( xcmsSet@phenoData[,1]==phenoDataClass[i])
    }

    N <- length(files)
    TIC <- vector("list",N)

    for (i in 1:N) {
        if (!is.null(xcmsSet) && rt == "corrected")
            rtcor <- xcmsSet@rt$corrected[[i]] else
        rtcor <- NULL
        TIC[[i]] <- getTIC(files[i],rtcor=rtcor)

        exportTicBpcTabular(TIC[[i]], files[i], "TIC", rt=rt)
    }

    pdf(pdfname,w=16,h=10)
    cols <- rainbow(N)
    lty = 1:N
    pch = 1:N
    #search for max x and max y in TICs
    xlim = range(sapply(TIC, function(x) range(x[,1])))
    ylim = range(sapply(TIC, function(x) range(x[,2])))
    ylim = c(-ylim[2], ylim[2])


    ##plot start
    if (length(phenoDataClass)>2){
        for (k in 1:(length(phenoDataClass)-1)){
            for (l in (k+1):length(phenoDataClass)){
                #print(paste(phenoDataClass[k],"vs",phenoDataClass[l],sep=" "))
                plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",phenoDataClass[k]," vs ",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
                colvect<-NULL
                for (j in 1:length(classnames[[k]])) {
                    tic <- TIC[[classnames[[k]][j]]]
                    # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
                    points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[k]][j]])
                }
                for (j in 1:length(classnames[[l]])) {
                    # i=class2names[j]
                    tic <- TIC[[classnames[[l]][j]]]
                    points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
                    colvect<-append(colvect,cols[classnames[[l]][j]])
                }
                legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)
            }
        }
    }#end if length >2
    if (length(phenoDataClass)==2){
        k=1
        l=2

        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",phenoDataClass[k],"vs",phenoDataClass[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
        colvect<-NULL
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }
        for (j in 1:length(classnames[[l]])) {
            # i=class2names[j]
            tic <- TIC[[classnames[[l]][j]]]
            points(tic[,1]/60, -tic[,2], col = cols[classnames[[l]][j]], pch = pch[classnames[[l]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[l]][j]])
        }
        legend("topright",paste(basename(files[c(classnames[[k]],classnames[[l]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==2

    #case where only one class
    if (length(phenoDataClass)==1){
        k=1
        ylim = range(sapply(TIC, function(x) range(x[,2])))

        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",phenoDataClass[k], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
        colvect<-NULL
        for (j in 1:length(classnames[[k]])) {
            tic <- TIC[[classnames[[k]][j]]]
            # points(tic[,1]/60, tic[,2], col = cols[i], pch = pch[i], type="l")
            points(tic[,1]/60, tic[,2], col = cols[classnames[[k]][j]], pch = pch[classnames[[k]][j]], type="l")
            colvect<-append(colvect,cols[classnames[[k]][j]])
        }

        legend("topright",paste(basename(files[c(classnames[[k]])])), col = colvect, lty = lty, pch = pch)

    }#end length ==1

    dev.off() #pdf(pdfname,w=16,h=10)

    invisible(TIC)
}



##
##  Get the polarities from all the samples of a condition
#@author Misharl Monsoor misharl.monsoor@sb-roscoff.fr ABiMS TEAM
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr ABiMS TEAM
getSampleMetadata <- function(xcmsSet=NULL, sampleMetadataOutput="sampleMetadata.tsv") {
    cat("Creating the sampleMetadata file...\n")

    #Create the sampleMetada dataframe
    sampleMetadata=xset@phenoData
    sampleNamesOrigin=rownames(sampleMetadata)
    sampleNamesMakeNames=make.names(sampleNamesOrigin)

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

    sampleMetadata$sampleMetadata=sampleNamesMakeNames
    sampleMetadata=cbind(sampleMetadata["sampleMetadata"],sampleMetadata["class"]) #Reorder columns
    rownames(sampleMetadata)=NULL

    #Create a list of files name in the current directory
    list_files=xset@filepaths
    #For each sample file, the following actions are done
    for (file in list_files){
        #Check if the file is in the CDF format
        if (!mzR:::netCDFIsFile(file)){

            # If the column isn't exist, with add one filled with NA
            if (is.null(sampleMetadata$polarity)) sampleMetadata$polarity=NA

            #Create a simple xcmsRaw object for each sample
            xcmsRaw=xcmsRaw(file)
            #Extract the polarity (a list of polarities)
            polarity=xcmsRaw@polarity
            #Verify if all the scans have the same polarity
            uniq_list=unique(polarity)
            if (length(uniq_list)>1){
                polarity="mixed"
            } else {
                polarity=as.character(uniq_list)
            }
            #Transforms the character to obtain only the sample name
            filename=basename(file)
            library(tools)
            samplename=file_path_sans_ext(filename)

            #Set the polarity attribute
            sampleMetadata$polarity[sampleMetadata$sampleMetadata==samplename]=polarity

            #Delete xcmsRaw object because it creates a bug for the fillpeaks step
            rm(xcmsRaw)
        }

    }

    write.table(sampleMetadata, sep="\t", quote=FALSE, row.names=FALSE, file=sampleMetadataOutput)

    return(list("sampleNamesOrigin"=sampleNamesOrigin,"sampleNamesMakeNames"=sampleNamesMakeNames))

}


##
## This function check if xcms will found all the files
##
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr ABiMS TEAM
checkFilesCompatibilityWithXcms <- function(directory) {
    cat("Checking files filenames compatibilities with xmcs...\n")
    # WHAT XCMS WILL FIND
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]","[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""),collapse = "|")
    info <- file.info(directory)
    listed <- list.files(directory[info$isdir], pattern = filepattern,recursive = TRUE, full.names = TRUE)
    files <- c(directory[!info$isdir], listed)
    files_abs <- file.path(getwd(), files)
    exists <- file.exists(files_abs)
    files[exists] <- files_abs[exists]
    files[exists] <- sub("//","/",files[exists])

    # WHAT IS ON THE FILESYSTEM
    filesystem_filepaths=system(paste("find $PWD/",directory," -not -name '\\.*' -not -path '*conda-env*' -type f -name \"*\"", sep=""), intern=T)
    filesystem_filepaths=filesystem_filepaths[grep(filepattern, filesystem_filepaths, perl=T)]

    # COMPARISON
    if (!is.na(table(filesystem_filepaths %in% files)["FALSE"])) {
        write("\n\nERROR: List of the files which will not be imported by xcmsSet",stderr())
        write(filesystem_filepaths[!(filesystem_filepaths %in% files)],stderr())
        stop("\n\nERROR: One or more of your files will not be import by xcmsSet. It may due to bad characters in their filenames.")
    }
}



##
## This function check if XML contains special caracters. It also checks integrity and completness.
##
#@author Misharl Monsoor misharl.monsoor@sb-roscoff.fr ABiMS TEAM
checkXmlStructure <- function (directory) {
    cat("Checking XML structure...\n")

    cmd=paste("IFS=$'\n'; for xml in $(find",directory,"-not -name '\\.*' -not -path '*conda-env*' -type f -iname '*.*ml*'); do if [ $(xmllint --nonet --noout \"$xml\" 2> /dev/null; echo $?) -gt 0 ]; then echo $xml;fi; done;")
    capture=system(cmd,intern=TRUE)

    if (length(capture)>0){
        #message=paste("The following mzXML or mzML file is incorrect, please check these files first:",capture)
        write("\n\nERROR: The following mzXML or mzML file(s) are incorrect, please check these files first:", stderr())
        write(capture, stderr())
        stop("ERROR: xcmsSet cannot continue with incorrect mzXML or mzML files")
    }

}


##
## This function check if XML contain special characters
##
#@author Misharl Monsoor misharl.monsoor@sb-roscoff.fr ABiMS TEAM
deleteXmlBadCharacters<- function (directory) {
    cat("Checking Non ASCII characters in the XML...\n")

    processed=F
    l=system( paste("find",directory, "-not -name '\\.*' -not -path '*conda-env*' -type f -iname '*.*ml*'"),intern=TRUE)
    for (i in l){
        cmd=paste("LC_ALL=C grep '[^ -~]' \"",i,"\"",sep="")
        capture=suppressWarnings(system(cmd,intern=TRUE))
        if (length(capture)>0){
            cmd=paste("perl -i -pe 's/[^[:ascii:]]//g;'",i)
            print( paste("WARNING: Non ASCII characters have been removed from the ",i,"file") )
            c=system(cmd,intern=TRUE)
            capture=""
            processed=T
        }
    }
    if (processed) cat("\n\n")
    return(processed)
}


##
## This function will compute MD5 checksum to check the data integrity
##
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getMd5sum <- function (directory) {
    cat("Compute md5 checksum...\n")
    # WHAT XCMS WILL FIND
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]","[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""),collapse = "|")
    info <- file.info(directory)
    listed <- list.files(directory[info$isdir], pattern = filepattern,recursive = TRUE, full.names = TRUE)
    files <- c(directory[!info$isdir], listed)
    exists <- file.exists(files)
    files <- files[exists]

    library(tools)

    #cat("\n\n")

    return(as.matrix(md5sum(files)))
}


# This function get the raw file path from the arguments
getRawfilePathFromArguments <- function(singlefile, zipfile, listArguments) {
    if (!is.null(listArguments[["zipfile"]]))           zipfile = listArguments[["zipfile"]]
    if (!is.null(listArguments[["zipfilePositive"]]))   zipfile = listArguments[["zipfilePositive"]]
    if (!is.null(listArguments[["zipfileNegative"]]))   zipfile = listArguments[["zipfileNegative"]]

    if (!is.null(listArguments[["singlefile_galaxyPath"]])) {
        singlefile_galaxyPaths = listArguments[["singlefile_galaxyPath"]];
        singlefile_sampleNames = listArguments[["singlefile_sampleName"]]
    }
    if (!is.null(listArguments[["singlefile_galaxyPathPositive"]])) {
        singlefile_galaxyPaths = listArguments[["singlefile_galaxyPathPositive"]];
        singlefile_sampleNames = listArguments[["singlefile_sampleNamePositive"]]
    }
    if (!is.null(listArguments[["singlefile_galaxyPathNegative"]])) {
        singlefile_galaxyPaths = listArguments[["singlefile_galaxyPathNegative"]];
        singlefile_sampleNames = listArguments[["singlefile_sampleNameNegative"]]
    }
    if (exists("singlefile_galaxyPaths")){
        singlefile_galaxyPaths = unlist(strsplit(singlefile_galaxyPaths,","))
        singlefile_sampleNames = unlist(strsplit(singlefile_sampleNames,","))

        singlefile=NULL
        for (singlefile_galaxyPath_i in seq(1:length(singlefile_galaxyPaths))) {
            singlefile_galaxyPath=singlefile_galaxyPaths[singlefile_galaxyPath_i]
            singlefile_sampleName=singlefile_sampleNames[singlefile_galaxyPath_i]
            singlefile[[singlefile_sampleName]] = singlefile_galaxyPath
        }
    }
    for (argument in c("zipfile","zipfilePositive","zipfileNegative","singlefile_galaxyPath","singlefile_sampleName","singlefile_galaxyPathPositive","singlefile_sampleNamePositive","singlefile_galaxyPathNegative","singlefile_sampleNameNegative")) {
        listArguments[[argument]]=NULL
    }
    return(list(zipfile=zipfile, singlefile=singlefile, listArguments=listArguments))
}


# This function retrieve the raw file in the working directory
#   - if zipfile: unzip the file with its directory tree
#   - if singlefiles: set symlink with the good filename
retrieveRawfileInTheWorkingDirectory <- function(singlefile, zipfile) {
    if(!is.null(singlefile) && (length("singlefile")>0)) {
        for (singlefile_sampleName in names(singlefile)) {
            singlefile_galaxyPath = singlefile[[singlefile_sampleName]]
            if(!file.exists(singlefile_galaxyPath)){
                error_message=paste("Cannot access the sample:",singlefile_sampleName,"located:",singlefile_galaxyPath,". Please, contact your administrator ... if you have one!")
                print(error_message); stop(error_message)
            }

            file.symlink(singlefile_galaxyPath,singlefile_sampleName)
        }
        directory = "."

    }
    if(!is.null(zipfile) && (zipfile!="")) {
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

    }
    return (directory)
}

# lib.r version="2.0.1"
#Authors ABiMS TEAM
#Lib.r for Galaxy Workflow4Metabo
#version 2.2
#Based on lib.r 2.1
#Modifications made by Guitton Yann 
#correct bug in Base Peak Chromatogram (BPC) option, not only TIC when scanrange used in xcmsSet
#Note if scanrange is used a warning is prompted in R console but do not stop PDF generation




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

  class<-as.vector(levels(xcmsSet@phenoData[,1])) #sometime phenoData have more than 1 column use first as class

  classnames<-vector("list",length(class))
  for (i in 1:length(class)){
    classnames[[i]]<-which( xcmsSet@phenoData[,1]==class[i])
  }

  N <- dim(phenoData(xcmsSet))[1]

  TIC <- vector("list",N)


  for (j in 1:N) {

    TIC[[j]] <- getBPC(files[j])
    #good for raw 
    # seems strange for corrected
    #errors if scanrange used in xcmsSetgeneration
    if (!is.null(xcmsSet) && rt == "corrected")
    rtcor <- xcmsSet@rt$corrected[[j]] else
    rtcor <- NULL
    
    TIC[[j]] <- getBPC(files[j],rtcor=rtcor)
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
  
  if (length(class)>2){
    for (k in 1:(length(class)-1)){
      for (l in (k+1):length(class)){
        #print(paste(class[k],"vs",class[l],sep=" ")) 
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k]," vs ",class[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")
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

  if (length(class)==2){
    k=1
    l=2
    colvect<-NULL
    plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Base Peak Chromatograms \n","BPCs_",class[k],"vs",class[l], sep=""), xlab = "Retention Time (min)", ylab = "BPC")

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

  class<-as.vector(levels(xcmsSet@phenoData[,1])) #sometime phenoData have more than 1 column use first as class

  classnames<-vector("list",length(class))
  for (i in 1:length(class)){
    classnames[[i]]<-which( xcmsSet@phenoData[,1]==class[i])
  }
  
  N <- length(files)
  TIC <- vector("list",N)

  for (i in 1:N) {
    if (!is.null(xcmsSet) && rt == "corrected")
      rtcor <- xcmsSet@rt$corrected[[i]] else
    rtcor <- NULL
    TIC[[i]] <- getTIC(files[i],rtcor=rtcor)
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
  if (length(class)>2){
    for (k in 1:(length(class)-1)){
      for (l in (k+1):length(class)){
        #print(paste(class[k],"vs",class[l],sep=" ")) 
        plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",class[k]," vs ",class[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
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
  if (length(class)==2){
    k=1
    l=2

    plot(0, 0, type="n", xlim = xlim/60, ylim = ylim, main = paste("Total Ion Chromatograms \n","TICs_",class[k],"vs",class[l], sep=""), xlab = "Retention Time (min)", ylab = "TIC")
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


#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 070115_NmrBucketing2galaxy_v1.R
## Marie Tremblay-Franco
## MetaboHUB: The French Infrastructure for Metabolomics and Fluxomics
## www.metabohub.fr/en
## marie.tremblay-franco@toulouse.inra.fr

runExampleL <- FALSE


##------------------------------
## Options
##------------------------------
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)


##------------------------------
## Libraries laoding
##------------------------------
# For parseCommandArgs function
library(batch) 

# R script call
source_local <- function(fname)
{
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
	source(paste(base_dir, fname, sep="/"))
}
#Import the different functions
source_local("NmrNormalization_script.R")

##------------------------------
## Errors ?????????????????????
##------------------------------


##------------------------------
## Constants
##------------------------------
topEnvC <- environment()
flagC <- "\n"


##------------------------------
## Script
##------------------------------
if(!runExampleL)
    argLs <- parseCommandArgs(evaluate=FALSE)


## Parameters Loading
##-------------------
  # Inputs
data <- read.table(argLs[["dataMatrix"]],check.names=FALSE,header=TRUE,sep="\t")
rownames(data) <- data[,1]
data <- data[,-1]

scaling <- argLs[["scalingMethod"]]
graphique <- argLs[["graphType"]]

if (scaling=='PQN')
{
 metadataSample <- read.table(argLs[["sampleMetadata"]],check.names=FALSE,header=TRUE,sep="\t")
 factor<- argLs[["factor"]]
 ControlGroup <- argLs[["controlGroup"]]
}
if (scaling=='QuantitativeVariable')
{
  metadataSample <- read.table(argLs[["sampleMetadata"]],check.names=FALSE,header=TRUE,sep="\t")
  factor <- argLs[["factor"]]
}

  # Outputs
nomGraphe <- argLs[["graphOut"]]
dataMatrixOut <- argLs[["dataMatrixOut"]]
sampleMetadataOut <- argLs[["sampleMetadataOut"]]
variableMetadataOut <- argLs[["variableMetadataOut"]]
log <- argLs[["logOut"]]

## Checking arguments
##-------------------
error.stock <- "\n"

if(length(error.stock) > 1)
  stop(error.stock)
  
  
## Computation
##------------
NormalizationResults <- NmrNormalization(dataMatrix=data,scalingMethod=scaling,sampleMetadata=metadataSample,
                                    bioFactor=factor,ControlGroup=ControlGroup,
                                    graph=graphique,nomFichier=nomGraphe,savLog.txtC=log)

data_normalized <- NormalizationResults[[1]]
data_sample <- NormalizationResults[[2]]
data_variable <- NormalizationResults[[3]]



## Saving
##-------
  # Data
data_normalized <- cbind(rownames(data_normalized),data_normalized)
colnames(data_normalized) <- c("Bucket",colnames(data_normalized)[-1])
write.table(data_normalized,file=argLs$dataMatrixOut,quote=FALSE,row.names=FALSE,sep="\t")
  # Sample
data_sample <- cbind(rownames(data_sample),data_sample)
colnames(data_sample) <- c("Sample",colnames(data_sample)[-1])
write.table(data_sample,file=argLs$sampleMetadataOut,quote=FALSE,row.names=FALSE,sep="\t")
  # Variable
data_variable <- cbind(rownames(data_variable),data_variable)
colnames(data_variable) <- c("Bucket",colnames(data_variable)[-1])
write.table(data_variable,file=argLs$variableMetadataOut,quote=FALSE,row.names=FALSE,sep="\t")


## Ending
##---------------------

cat("\nEnd of 'NMR Normalization' Galaxy module call: ", as.character(Sys.time()), sep = "")

## sink(NULL)

options(stringsAsFactors = strAsFacL)

rm(list = ls())

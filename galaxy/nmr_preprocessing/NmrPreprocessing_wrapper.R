#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 170116_NmrPreprocessing.R
## Manon Martin and Marie Tremblay-Franco

##======================================================
##======================================================
# Preamble
##======================================================
##======================================================

runExampleL <- FALSE

##------------------------------
## Options
##------------------------------
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

##------------------------------
## Libraries laoding
##------------------------------
<<<<<<< HEAD
library(batch)
library(ptw)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
=======
# For parseCommandArgs function
ip <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ip) <- NULL
ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
# print(ip, row.names=FALSE)

# library(batch)
# library(ptw)
# library(Matrix)
# library(ggplot2)
# library(gridExtra)
# library(reshape2)
>>>>>>> 3cdb294d8554f70467930c415d2a905d833fbcec


# R script call
source_local <- function(fname)
{
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
	source(paste(base_dir, fname, sep="/"))
}
#Import the different functions
source_local("NmrPreprocessing_script.R")
source_local("DrawFunctions.R")

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

# log file
print(argLs[["logOut"]])

## Starting
cat("\nStart of 'Preprocessing' Galaxy module call: ", as.character(Sys.time()), sep = "")


##======================================================
##======================================================
## Parameters Loading
##======================================================
##======================================================

# graphical inputs
FirstOPCGraph <- argLs[["FirstOPCGraph"]]
SSGraph <- argLs[["SSGraph"]]
ApodGraph <- argLs[["ApodGraph"]]
FTGraph <- argLs[["FTGraph"]]
SRGraph <- argLs[["SRGraph"]]
ZeroOPCGraph <- argLs[["ZeroOPCGraph"]]
BCGraph <- argLs[["BCGraph"]]
FinalGraph <- argLs[["FinalGraph"]]


# 1rst order phase correction ------------------------
  # Inputs
	## Data matrix
Fid_data0 <- read.table(argLs[["dataMatrixFid"]],header=TRUE, check.names=FALSE, sep='\t')
# Fid_data0 <- Fid_data0[,-1]
Fid_data0 <- as.matrix(Fid_data0)

	## Samplemetadata
samplemetadataFid <- read.table(argLs[["sampleMetadataFid"]],check.names=FALSE,header=TRUE,sep="\t")
samplemetadataFid <- as.matrix(samplemetadataFid)


# water and solvent(s) correction ------------------------
  # Inputs
lambda <- argLs[["lambda"]]
ptwSS1 <- argLs[["ptwSS"]]
ptwSS <- FALSE
if (ptwSS1=="YES")
 	ptwSS <- TRUE
<<<<<<< HEAD

=======
	
>>>>>>> 3cdb294d8554f70467930c415d2a905d833fbcec

# apodization -----------------------------------------
  # Inputs
phase=0
rectRatio=1/2
gaussLB=1
expLB=1
apodization <- argLs[["apodizationMethod"]]

if (apodization=='exp'){
  expLB <- argLs[["expLB"]]
  } else if (apodization=='cos2'){
  phase <- argLs[["phase"]]
  } else if (apodization=='hanning'){
  phase <- argLs[["phase"]]
  } else if (apodization=='hamming'){
  phase <- argLs[["phase"]]
  } else if (apodization=='blockexp'){
  rectRatio <- argLs[["rectRatio"]]
  expLB <- argLs[["expLB"]]
  } else if (apodization=='blockcos2'){
  rectRatio <- argLs[["rectRatio"]]
  } else if (apodization=='gauss'){
  rectRatio <- argLs[["rectRatio"]]
  gaussLB <- argLs[["gaussLB"]]
  }		


# Fourier transform ----------------------------------
  # Inputs


# Internal referencering ----------------------------------
  # Inputs
<<<<<<< HEAD
shiftTreshold = 2 # c
ppm = TRUE
shiftReferencingRangeList = NULL  # fromto.RC
pctNear0 = 0.02 # pc 
rowindex_graph = NULL
ppm_ref = 0 # ppm.ref

# 
# shiftReferencing <- argLs[["shiftReferencing"]]
# print(shiftReferencing)
# 
# if (shiftReferencing=="YES")
# {
#   
	# shiftReferencingMethod <- argLs[["shiftReferencingMethod"]]
	# 
	# if (shiftReferencingMethod == "thres")	{
	# 	shiftTreshold <- argLs[["shiftTreshold"]]
	# }
	
	shiftReferencingRange <- argLs[["shiftReferencingRange"]]
	
	if (shiftReferencingRange == "near0"){
	  pctNear0 <- argLs[["pctNear0"]]
	}
	
	if (shiftReferencingRange == "window"){
	  shiftReferencingRangeList <- list()
	  for(i in which(names(argLs)=="shiftReferencingRangeLeft")) 
	  {
  		shiftReferencingRangeLeft <- argLs[[i]]
  		shiftReferencingRangeRight <- argLs[[i+1]]
  		shiftReferencingRangeList <- c(shiftReferencingRangeList,list(c(shiftReferencingRangeLeft,shiftReferencingRangeRight)))
	  }
	}
	
	shiftHandling <- argLs[["shiftHandling"]]
	
=======
ShiftReferencing <- argLs[["shift_referencing.choice"]]
print(ShiftReferencing)
if (ShiftReferencing=="YES")
{
	shiftTreshold = 2 # c
	ppmxaxis = TRUE
	shiftReferencingRangeList = NULL  # fromto.TMSP
	pctNear0 = 0.02 # pc 
	rowindex_graph = NULL
	shiftReferencingMethod <- argLs[["shiftReferencingMethod"]]
	if (shiftReferencingMethod == "thres")
	{
		shiftTreshold <- argLs[["shiftTreshold"]]
	}
	shiftReferencingRange <- argLs[["shiftReferencingRange"]]
	print(shiftReferencingRange)
	if (shiftReferencingRange == "near0"){
		pctNear0 <- argLs[["pctNear0"]]
	}
	 if (shiftReferencingRange == "window"){
			shiftReferencingRangeList <- NULL
			shiftReferencingRangeLeft <- argLs[["shiftReferencingRangeLeft"]]
			shiftReferencingRangeRight <- argLs[["shiftReferencingRangeRight"]]
			shiftReferencingRangeList <- list(shiftReferencingRangeList,c(shiftReferencingRangeLeft,shiftReferencingRangeRight))
	 }
	shiftHandling <- argLs[["shiftHandling"]]
	shiftReferencing <- argLs[["shiftReferencing"]]
	if (shiftReferencingRange == "near0"){
		  pctNear0 <- argLs[["pctNear0"]]
	}
		
	if (shiftReferencingRange == "window"){
		shiftReferencingRangeList <- NULL
		shiftReferencingRangeLeft <- argLs[["shiftReferencingRangeLeft"]]
		shiftReferencingRangeRight <- argLs[["shiftReferencingRangeRight"]]
		shiftReferencingRangeList <- list(shiftReferencingRangeList,c(shiftReferencingRangeLeft,shiftReferencingRangeRight))
	}
	shiftHandling <- argLs[["shiftHandling"]]
}
>>>>>>> 3cdb294d8554f70467930c415d2a905d833fbcec

	
# }


# Zero Order Phase Correction -------------------------------
  # Inputs
plot_rms = NULL
returnAngle = FALSE
createWindow = TRUE
angle = NULL
plot_spectra = FALSE
ppm = TRUE
exclude = NULL
zeroOrderPhaseMethod <- argLs[["zeroOrderPhaseMethod"]]
										   
if (zeroOrderPhaseMethod=='manual'){
  angle <- argLs[["angle"]]
}

excludeZoneZeroPhase <- argLs[["excludeZoneZeroPhase.choice"]]
if (excludeZoneZeroPhase == 'YES') {
  excludeZoneZeroPhaseList <- list()
  for(i in which(names(argLs)=="excludeZoneZeroPhase_left")) {
    excludeZoneZeroPhaseLeft <- argLs[[i]]
    excludeZoneZeroPhaseRight <- argLs[[i+1]]
    excludeZoneZeroPhaseList <- c(excludeZoneZeroPhaseList,list(c(excludeZoneZeroPhaseLeft,excludeZoneZeroPhaseRight)))
  }
  exclude <- excludeZoneZeroPhaseList
}

# Baseline Correction -------------------------------
  # Inputs
ptwBc <- as.logical(argLs[["ptwBc"]])
maxIter <- argLs[["maxIter"]] 
lambdaBc <- argLs[["lambdaBc"]] 
pBc <- argLs[["pBc"]] 
epsilon <- argLs[["epsilon"]] 
<<<<<<< HEAD

=======
>>>>>>> 3cdb294d8554f70467930c415d2a905d833fbcec


# transformation of negative values -------------------------------
  # Inputs
NegativetoZero <- argLs[["NegativetoZero"]]


  # Outputs
nomGraphe <- argLs[["graphOut"]]
dataMatrixOut <- argLs[["dataMatrixOut"]]
log <- argLs[["logOut"]]



## Checking arguments
##-------------------
error.stock <- "\n"

if(length(error.stock) > 1)
  stop(error.stock)
  

##======================================================
##======================================================  
## Computation
##======================================================
##======================================================

pdf(nomGraphe, onefile = TRUE, width = 13, height = 13)

# FirstOrderPhaseCorrection ---------------------------------
Fid_data <- FirstOrderPhaseCorrection(Fid_data0, Fid_info = samplemetadataFid, group_delay = NULL)

if (FirstOPCGraph == "YES") {
  title = "FIDs after First Order Phase Correction"
  DrawSignal(Fid_data, subtype = "stacked",
             ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
             xlab = "Frequency", num.stacked = 4, 
             main = title, createWindow=FALSE)
}

# SolventSuppression ---------------------------------
Fid_data <- SolventSuppression(Fid_data, lambda.ss = lambda, ptw.ss = ptwSS, plotSolvent = F, returnSolvent = F)
	
if (SSGraph == "YES") {
  title = "FIDs after Solvent Suppression "
  DrawSignal(Fid_data, subtype = "stacked",
             ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
             xlab = "Frequency", num.stacked = 4, 
             main = title, createWindow=FALSE)
}


# Apodization ---------------------------------	
Fid_data <- Apodization(Fid_data, Fid_info = samplemetadataFid, DT = NULL, 
                         type.apod = apodization, phase = phase, rectRatio = rectRatio, gaussLB = gaussLB, expLB = expLB, plotWindow = F, returnFactor = F)

if (ApodGraph == "YES") {
  title = "FIDs after Apodization"
  DrawSignal(Fid_data, subtype = "stacked",
             ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
             xlab = "Frequency", num.stacked = 4, 
             main = title, createWindow=FALSE)
}


# FourierTransform ---------------------------------
Spectrum_data <- FourierTransform(Fid_data, Fid_info = samplemetadataFid, reverse.axis = TRUE)

<<<<<<< HEAD

=======
>>>>>>> 3cdb294d8554f70467930c415d2a905d833fbcec
if (FTGraph == "YES") {
  title = "Fourier transformed spectra"
  DrawSignal(Spectrum_data, subtype = "stacked",
             ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
             xlab = "Frequency", num.stacked = 4, 
             main = title, createWindow=FALSE)
}


# InternalReferencing ---------------------------------
<<<<<<< HEAD
# if (shiftReferencing=="YES") {
Spectrum_data <- InternalReferencing(Spectrum_data, samplemetadataFid, method = "max", range = shiftReferencingRange,
                                     ppm.ref = 0, shiftHandling = shiftHandling,ppm = TRUE,
									 c = shiftTreshold, fromto.RC = shiftReferencingRangeList, pc = pctNear0)

if (SRGraph == "YES") {
  title = "Spectra after Shift Referencing"
  DrawSignal(Spectrum_data, subtype = "stacked",
             ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
             xlab = "Frequency", num.stacked = 4, 
             main = title, createWindow=FALSE)
=======
if (ShiftReferencing=="YES")
{
	Spectrum_data <- InternalReferencing(Spectrum_data, samplemetadataFid, method = shiftReferencingMethod, range = shiftReferencingRange, shiftHandling = shiftHandling,
										c = shiftTreshold, fromto.TMSP = shiftReferencingRangeList, pc = pctNear0)
	
	title = "Spectra after Shift referencing"
	DrawSignal(Spectrum_data, subtype = "stacked",
               ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
               xlab = "Frequency", num.stacked = 4, 
               main.title = title, createWindow=FALSE)
>>>>>>> 3cdb294d8554f70467930c415d2a905d833fbcec
}


# ZeroOrderPhaseCorrection ---------------------------------
Spectrum_data  <- ZeroOrderPhaseCorrection(Spectrum_data, method = zeroOrderPhaseMethod,
                                           plot_rms = plot_rms, returnAngle = returnAngle,
                                           createWindow = createWindow,angle = angle,
                                           plot_spectra = plot_spectra,
                                           ppm = ppm, exclude = exclude)

if (ZeroOPCGraph == "YES") {
title = "Spectra after Zero Order Phase Correction"
DrawSignal(Spectrum_data, subtype = "stacked",
           ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
           xlab = "Frequency", num.stacked = 4, 
<<<<<<< HEAD
           main = title, createWindow=FALSE)
}

# BaselineCorrection ---------------------------------									 
Spectrum_data <- BaselineCorrection(Spectrum_data, ptw.bc = ptwBc, maxIter = maxIter, lambda.bc = lambdaBc, p.bc = pBc, eps = epsilon, returnBaseline = F) 

if (BCGraph == "YES") {
title = "Spectra after Baseline Correction"
DrawSignal(Spectrum_data, subtype = "stacked",
           ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
           xlab = "Frequency", num.stacked = 4, 
           main = title, createWindow=FALSE)
}

=======
           main.title = title, createWindow=FALSE)


# BaselineCorrection ---------------------------------			
Spectrum_data <- BaselineCorrection(Spectrum_data, ptw.bc = ptwBc, maxIter = maxIter, lambda.bc = lambdaBc, p.bc = pBc, eps = epsilon, returnBaseline = F) 

title = "Spectra after Baseline Correction"
DrawSignal(Spectrum_data, subtype = "stacked",
			ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
			xlab = "Frequency", num.stacked = 4, 
			main.title = title, createWindow=FALSE)
 	
>>>>>>> 3cdb294d8554f70467930c415d2a905d833fbcec

# NegativeValuesZeroing ---------------------------------
if (NegativetoZero=="YES") {
  Spectrum_data <- NegativeValuesZeroing(Spectrum_data)
}

<<<<<<< HEAD
if (FinalGraph == "YES") {
  title = "Final preprocessed spectra"
  DrawSignal(Spectrum_data, subtype = "stacked",
             ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T, 
             xlab = "Frequency", num.stacked = 4, 
             main = title, createWindow=FALSE)
}

=======
>>>>>>> 3cdb294d8554f70467930c415d2a905d833fbcec
invisible(dev.off())


##======================================================
##======================================================
## Saving
##======================================================
##======================================================
# Data Matrix
<<<<<<< HEAD
write.table(t(Re(Spectrum_data)),file=argLs$dataMatrixOut, quote=FALSE, row.names=TRUE, sep="\t", col.names=TRUE)
=======
write.table(Spectrum_data,file=argLs$dataMatrixOut, quote=FALSE, row.names=TRUE, sep="\t", col.names=TRUE)
>>>>>>> 3cdb294d8554f70467930c415d2a905d833fbcec


## Ending
cat("\nEnd of 'Preprocessing' Galaxy module call: ", as.character(Sys.time()), sep = "")
options(stringsAsFactors = strAsFacL)
rm(list = ls())

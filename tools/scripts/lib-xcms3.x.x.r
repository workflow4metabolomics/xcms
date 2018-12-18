

#@TODO: remove this function as soon as we can use xcms 3.x.x from Bioconductor 3.7
# https://github.com/sneumann/xcms/issues/250
groupnamesW4M <- function(xdata, mzdec = 0, rtdec = 0) {
    mzfmt <- paste("%.", mzdec, "f", sep = "")
    rtfmt <- paste("%.", rtdec, "f", sep = "")

    gnames <- paste("M", sprintf(mzfmt, featureDefinitions(xdata)[,"mzmed"]), "T",
                    sprintf(rtfmt, featureDefinitions(xdata)[,"rtmed"]), sep = "")

    if (any(dup <- duplicated(gnames)))
        for (dupname in unique(gnames[dup])) {
            dupidx <- which(gnames == dupname)
            gnames[dupidx] <- paste(gnames[dupidx], seq(along = dupidx), sep = "_")
        }

    return (gnames)
}

#@TODO: remove this function as soon as we can use xcms 3.x.x from Bioconductor 3.7
# https://github.com/sneumann/xcms/issues/247
.concatenate_XCMSnExp <- function(...) {
    x <- list(...)
    if (length(x) == 0)
        return(NULL)
    if (length(x) == 1)
        return(x[[1]])
    ## Check that all are XCMSnExp objects.
    if (!all(unlist(lapply(x, function(z) is(z, "XCMSnExp")))))
        stop("All passed objects should be 'XCMSnExp' objects")
    new_x <- as(.concatenate_OnDiskMSnExp(...), "XCMSnExp")
    ## If any of the XCMSnExp has alignment results or detected features drop
    ## them!
    x <- lapply(x, function(z) {
        if (hasAdjustedRtime(z)) {
            z <- dropAdjustedRtime(z)
            warning("Adjusted retention times found, had to drop them.")
        }
        if (hasFeatures(z)) {
            z <- dropFeatureDefinitions(z)
            warning("Feature definitions found, had to drop them.")
        }
        z
    })
    ## Combine peaks
    fls <- lapply(x, fileNames)
    startidx <- cumsum(lengths(fls))
    pks <- lapply(x, chromPeaks)
    procH <- lapply(x, processHistory)
    for (i in 2:length(fls)) {
        pks[[i]][, "sample"] <- pks[[i]][, "sample"] + startidx[i - 1]
        procH[[i]] <- lapply(procH[[i]], function(z) {
            z@fileIndex <- as.integer(z@fileIndex + startidx[i - 1])
            z
            })
    }
    pks <- do.call(rbind, pks)
    new_x@.processHistory <- unlist(procH)
    chromPeaks(new_x) <- pks
    if (validObject(new_x))
        new_x
}

#@TODO: remove this function as soon as we can use xcms 3.x.x from Bioconductor 3.7
# https://github.com/sneumann/xcms/issues/247
.concatenate_OnDiskMSnExp <- function(...) {
    x <- list(...)
    if (length(x) == 0)
        return(NULL)
    if (length(x) == 1)
        return(x[[1]])
    ## Check that all are XCMSnExp objects.
    if (!all(unlist(lapply(x, function(z) is(z, "OnDiskMSnExp")))))
        stop("All passed objects should be 'OnDiskMSnExp' objects")
    ## Check processingQueue
    procQ <- lapply(x, function(z) z@spectraProcessingQueue)
    new_procQ <- procQ[[1]]
    is_ok <- unlist(lapply(procQ, function(z)
        !is.character(all.equal(new_procQ, z))
        ))
    if (any(!is_ok)) {
        warning("Processing queues from the submitted objects differ! ",
                "Dropping the processing queue.")
        new_procQ <- list()
    }
    ## processingData
    fls <- lapply(x, function(z) z@processingData@files)
    startidx <- cumsum(lengths(fls))
    ## featureData
    featd <- lapply(x, fData)
    ## Have to update the file index and the spectrum names.
    for (i in 2:length(featd)) {
        featd[[i]]$fileIdx <- featd[[i]]$fileIdx + startidx[i - 1]
        rownames(featd[[i]]) <- MSnbase:::formatFileSpectrumNames(
                                              fileIds = featd[[i]]$fileIdx,
                                              spectrumIds = featd[[i]]$spIdx,
                                              nSpectra = nrow(featd[[i]]),
                                              nFiles = length(unlist(fls))
                                          )
    }
    featd <- do.call(rbind, featd)
    featd$spectrum <- 1:nrow(featd)
    ## experimentData
    expdata <- lapply(x, function(z) {
        ed <- z@experimentData
        data.frame(instrumentManufacturer = ed@instrumentManufacturer,
                   instrumentModel = ed@instrumentModel,
                   ionSource = ed@ionSource,
                   analyser = ed@analyser,
                   detectorType = ed@detectorType,
                   stringsAsFactors = FALSE)
    })
    expdata <- do.call(rbind, expdata)
    expdata <- new("MIAPE",
                   instrumentManufacturer = expdata$instrumentManufacturer,
                   instrumentModel = expdata$instrumentModel,
                   ionSource = expdata$ionSource,
                   analyser = expdata$analyser,
                   detectorType = expdata$detectorType)

    ## protocolData
    protodata <- lapply(x, function(z) z@protocolData)
    if (any(unlist(lapply(protodata, nrow)) > 0))
        warning("Found non-empty protocol data, but merging protocol data is",
                " currently not supported. Skipped.")
    ## phenoData
    pdata <- do.call(rbind, lapply(x, pData))
    res <- new(
        "OnDiskMSnExp",
        phenoData = new("NAnnotatedDataFrame", data = pdata),
        featureData = new("AnnotatedDataFrame", featd),
        processingData = new("MSnProcess",
                             processing = paste0("Concatenated [", date(), "]"),
                             files = unlist(fls), smoothed = NA),
        experimentData = expdata,
        spectraProcessingQueue = new_procQ)
    if (validObject(res))
        res
}

#@TODO: remove this function as soon as we can use xcms 3.x.x from Bioconductor 3.7
# https://github.com/sneumann/xcms/issues/247
c.XCMSnExp <- function(...) {
    .concatenate_XCMSnExp(...)
}

#@TODO: remove this function as soon as we can use xcms 3.x.x from Bioconductor 3.7
# https://github.com/sneumann/xcms/issues/247
c.MSnbase <- function(...) {
    .concatenate_OnDiskMSnExp(...)
}

#@TODO: remove this function as soon as we can use xcms 3.x.x from Bioconductor 3.7
# https://github.com/workflow4metabolomics/xcms/issues/117
.getChromPeakData <- function(object, peakArea, sample_idx,
                              mzCenterFun = "weighted.mean",
                              cn = c("mz", "rt", "into", "maxo", "sample")) {
    if (length(fileNames(object)) != 1)
        stop("'object' should be an XCMSnExp for a single file!")
    ncols <- length(cn)
    res <- matrix(ncol = ncols, nrow = nrow(peakArea))
    colnames(res) <- cn
    res[, "sample"] <- sample_idx
    res[, c("mzmin", "mzmax")] <- peakArea[, c("mzmin", "mzmax")]
    ## Load the data
    message("Requesting ", nrow(res), " peaks from ",
             basename(fileNames(object)), " ... ", appendLF = FALSE)
    spctr <- spectra(object, BPPARAM = SerialParam())
    mzs <- lapply(spctr, mz)
    valsPerSpect <- lengths(mzs)
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
    rm(spctr)
    mzs <- unlist(mzs, use.names = FALSE)
    mzs_range <- range(mzs)
    rtim <- rtime(object)
    rtim_range <- range(rtim)
    for (i in seq_len(nrow(res))) {
        rtr <- peakArea[i, c("rtmin", "rtmax")]
        mzr <- peakArea[i, c("mzmin", "mzmax")]
        ## If the rt range is completely out; additional fix for #267
        if (rtr[2] < rtim_range[1] | rtr[1] > rtim_range[2] |
            mzr[2] < mzs_range[1] | mzr[1] > mzs_range[2]) {
            res[i, ] <- rep(NA_real_, ncols)
            next
        }
        ## Ensure that the mz and rt region is within the range of the data.
        rtr[1] <- max(rtr[1], rtim_range[1])
        rtr[2] <- min(rtr[2], rtim_range[2])
        mzr[1] <- max(mzr[1], mzs_range[1])
        mzr[2] <- min(mzr[2], mzs_range[2])
        mtx <- .rawMat(mz = mzs, int = ints, scantime = rtim,
                       valsPerSpect = valsPerSpect, rtrange = rtr,
                       mzrange = mzr)
        if (length(mtx)) {
            if (any(!is.na(mtx[, 3]))) {
                ## How to calculate the area: (1)sum of all intensities / (2)by
                ## the number of data points (REAL ones, considering also NAs)
                ## and multiplied with the (3)rt width.
                ## (1) sum(mtx[, 3], na.rm = TRUE)
                ## (2) sum(rtim >= rtr[1] & rtim <= rtr[2]) - 1 ; if we used
                ## nrow(mtx) here, which would correspond to the non-NA
                ## intensities within the rt range we don't get the same results
                ## as e.g. centWave. Using max(1, ... to avoid getting Inf in
                ## case the signal is based on a single data point.
                ## (3) rtr[2] - rtr[1]
                res[i, "into"] <- sum(mtx[, 3], na.rm = TRUE) *
                    ((rtr[2] - rtr[1]) /
                     max(1, (sum(rtim >= rtr[1] & rtim <= rtr[2]) - 1)))
                maxi <- which.max(mtx[, 3])
                res[i, c("rt", "maxo")] <- mtx[maxi[1], c(1, 3)]
                res[i, c("rtmin", "rtmax")] <- rtr
                ## Calculate the intensity weighted mean mz
                meanMz <- do.call(mzCenterFun, list(mtx[, 2], mtx[, 3]))
                if (is.na(meanMz)) meanMz <- mtx[maxi[1], 2]
                res[i, "mz"] <- meanMz
            } else {
                res[i, ] <- rep(NA_real_, ncols)
            }
        } else {
            res[i, ] <- rep(NA_real_, ncols)
        }
    }
    message("got ", sum(!is.na(res[, "into"])), ".")
    return(res)
}

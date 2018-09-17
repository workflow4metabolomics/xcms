

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

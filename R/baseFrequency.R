# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  DNA Frequency Matrix Generator                  ><>< #
# ><>< ================================================================ ><>< #

baseCheck <- function(basename){
    basetag <- tolower(basename)
    invalid <- !(basetag %in% c("aa", "dna", "codon"))
    msg <- paste("Incorrect `basename` specification.",
        "\n\t Input one of `dna`, `codon` or `aa`.", sep="")
    if(invalid) stop(msg, call.=FALSE)
    return(basetag)
}

baseFrequency.noncodon <- function(baseseqs){
    countmatrix <- consensusMatrix(baseseqs, baseOnly=TRUE)
    excessrow <- which( rownames(countmatrix) == "other" )
    finalFrame <- countmatrix[-excessrow,]
    return(finalFrame)
}

codonCensusVector <- function(codonvector){
    outMatrix <- matrix(0,1,62)
    colnames(outMatrix) <- c(senseCodon, "other")
    chooSense <- codonvector %in% senseCodon
    codonvector[!chooSense] <- "other"
    codonTab <- table(codonvector)
    outMatrix[1,names(codonTab)] <- c(codonTab)
    return(outMatrix)
}

baseFrequency.codon <- function(bstringset){
    codonsequences <- bstringCodons(bstringset)
    finalFrame <- apply(codonsequences, 2, codonCensusVector)[-62,]
    rownames(finalFrame) <- senseCodon
    return(finalFrame)
}

baseScreen <- function(nleaf, freqMatrix, baseseqs){
    aptLength <- min(.5 * nrow(freqMatrix), 0.25 * max(width(baseseqs)))
    exitest <- sum(rowSums(freqMatrix) != 0) < aptLength
    exitbin <- sum(colSums(freqMatrix) < .9 * nleaf) > (.5 * ncol(freqMatrix))
    msg <- paste("Significant data imbalance detected.\nYou may ",
        "need to check that correct `basename` specified.",sep="")
    if(exitest | exitbin) warning(msg, call.=FALSE)
    return(0)
}

baseFrequency <- function(fastafile, basename){
    cleanbase <- baseCheck(basename)
    baseseqs <- switch(cleanbase,
        aa = readAAStringSet(fastafile),
        dna = readDNAStringSet(fastafile),
        codon = readBStringSet(fastafile)
    )
    freqMatrix <- switch(cleanbase,
        codon = baseFrequency.codon(baseseqs),
        aa = baseFrequency.noncodon(baseseqs),
        dna = baseFrequency.noncodon(baseseqs)
    )
    nleaf <- length(baseseqs)
    warns <- baseScreen(nleaf, freqMatrix, baseseqs)
    bType <- c(dna="DNA", codon="Codon", aa="Amino Acid")[cleanbase]
    bfoutput <- new("baseSummary", baseType=bType,
        extantSize=nleaf, seqsNumerics=freqMatrix)
    return(bfoutput)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                 Nucleotide to Codon Frequencies.                 ><>< #
# ><>< ================================================================ ><>< #

codon2int <- function(codon){
    nucIDs <- c(A=1, C=2, G=3, T=4)
    codonxters <- unlist( strsplit(codon, ""))
    nuConvert <- as.integer(nucIDs[codonxters])
    names(nuConvert) <- NULL
    return(nuConvert)
}

vetFreqs <- function(nucFracs){
    nTags <- c("A","C","G","T")
    init <- "\nNucleotide frequency vector "
    note3 <- paste(init, "must be named ",
        "properly or not contain `NA`.", sep="")
    note2 <- paste(init, "should only contain ",
        "non-negative numerical values.", sep="")
    note0 <- paste(init, "must have at least 4 elements.", sep="")
    note4 <- "\nOnly first four values in nucleotide frequency vector used."
    crime1 <- (length( names(nucFracs)) == 0) & (length(nucFracs) < 4)
    if( !all(nucFracs >= 0)) stop(note2, call.=FALSE)
    if(crime1) stop(note0, call.=FALSE)
    if( !is(names(nucFracs),"NULL")) nucFracs <- nucFracs[nTags]
    if( any( is.na(nucFracs))) stop(note3, call.=FALSE)
    if( length(nucFracs)>4) warning(note4, call.=FALSE)
    nucEdit <- nucFracs[seq(1,4)]
    newNucs <- nucEdit / sum(nucEdit)
    names(newNucs) <- nTags
    return(newNucs)
}

n2cFreqs <- function(nucFracs){
    nucFreqs <- vetFreqs(nucFracs)
    allcodonID <- vapply(senseCodon, codon2int, FUN.VALUE=vector("integer",3))
    codonWeights <- apply(allcodonID, 2, function(j) prod(nucFreqs[j]))
    codonFreqs <- codonWeights / sum(codonWeights)
    return(codonFreqs)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

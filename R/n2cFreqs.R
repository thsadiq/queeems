# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                 Nucleotide to Codon Frequencies.                 ><>< #
# ><>< ================================================================ ><>< #

codon2int <- function(codon){
    nucIDs <- c(A=1, C=2, G=3, T=4)
    codonxters <- unlist( strsplit(codon, ""))
    nuConvert <- nucIDs[codonxters]
    names(nuConvert) <- NULL
    return(nuConvert)
}

vetFreqs <- function(nucFracs){
    nTags <- c("A","C","G","T")
    init <- "\nNucleotide frequency vector "
    note0 <- paste0(init, "must have at least 4 elements.")
    note3 <- paste0(init, "must be named properly or not contain `NA`.")
    note2 <- paste0(init,"should only contain non-negative numerical values.")
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
    allcodonID <- sapply(senseCodon, codon2int)
    codonWeights <- apply(allcodonID, 2, function(j) prod(nucFreqs[j]))
    codonFreqs <- codonWeights / sum(codonWeights)
    return(codonFreqs)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

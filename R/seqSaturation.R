# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                        Saturation Indices                        ><>< #
# ><>< ================================================================ ><>< #

nucPIgen <- function(wsize, wless1, weightmax){
    allweights <- fubarweights(wsize, wless1, weightmax, FALSE)
    nucIndices <- permutations(wsize, 4)
    nucWeights <- apply(nucIndices, 1, function(i) allweights[i])
    nucFreqs <- apply(nucWeights, 2, softmax, maxadj=TRUE)
    return(nucFreqs)
}

siteLLmax <- function(basevector, pimatrix){
    loglikely <- apply(pimatrix, 2, function(j) mnomLogl(basevector,j))
    maxlikely <- 0
    if(!all(loglikely == 0)){
        postprob <- softmax(loglikely, TRUE)
        maxlikely <- max(postprob)
    }
    return(maxlikely)
}

seqsat1 <- function(fastafile, wsize, wless1, weightmax, bfcutoff, sspi){
    codoncounts <- basecensus( baseFrequency(fastafile, "codon"))
    nucleotidecounts <- basecensus( baseFrequency(fastafile, "dna"))
    nucProps <- nucPIgen(wsize, wless1, weightmax)
    codonProps <- apply(nucProps, 2, n2cFreqs)
    codonLL <- apply(codoncounts, 2, siteLLmax, pimatrix=codonProps)
    nucLLmax <- apply(nucleotidecounts, 2, siteLLmax, pimatrix=nucProps)
    meaNucID <- vapply(seq(1,3), function(i) seq(i,length(nucLLmax),by=3),
        FUN.VALUE=vector("numeric",length(nucLLmax)/3))
    nucLL <- apply(meaNucID, 1, function(i) sum(nucLLmax[i]))
    siteBF <- nucLL / codonLL
    checkBF <- (codonLL == 0) | (nucLL == 0)
    siteBF[checkBF] <- 0
    red <- which( sum( colSums(codoncounts) > 0) == 1)
    a1out <- new("saturateBF", bf.pv=siteBF, nuLL=codonLL, altLL=nucLL,
        noInfo=red, nsites=ncol(codoncounts), bfcutoff=bfcutoff, sspi=sspi)
    return( a1out )
}

scanapp <- function(approach){
    meth <- toupper(approach)
    wrd <- "\nIncorrect `approach` specified. Expects \U2208 {`A`,`B`,`C`}"
    if(!(meth %in% c("A","B","C"))) stop(wrd, call.=FALSE)
    return(meth)
}

seqSaturation <- function(fastafile, approach="A", wsize=6,
    wless1=0.40, weightmax=7, bfcutoff=3, sspi=0.01){
    method <- scanapp(approach)
    sanalyses <- switch(method,
        A = seqsat1(fastafile, wsize, wless1, weightmax, bfcutoff, sspi),
        B = stop("\nApproach B is yet to be finalised."),
        C = stop("\nApproach C is yet to be finalised."))
    return(sanalyses)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

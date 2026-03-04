# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                        Saturation Indices                        ><>< #
# ><>< ================================================================ ><>< #

nucPIgen <- function(wsize, wless1, weightmax){
    allweights <- fubarweights(wsize, wless1, weightmax, FALSE)
    nucIndices <- gtools::permutations(wsize, 4)
    nucWeights <- apply(nucIndices, 1, function(i) allweights[i])
    nucFreqs <- apply(nucWeights, 2, softmax, maxadj=TRUE)
    return(nucFreqs)
}

siteLLmax <- function(basevector, pimatrix){
    loglikely <- apply(pimatrix, 2, function(j) mnomLogl(basevector,j))
    maxlikely <- max(loglikely)
    return(maxlikely)
}

seqsat1 <- function(fastafile, wsize, wless1, weightmax){
    codoncounts <- basecensus( baseFrequency(fastafile, "codon"))
    nucleotidecounts <- basecensus( baseFrequency(fastafile, "dna"))
    nucProps <- nucPIgen(wsize, wless1, weightmax)
    codonProps <- apply(nucProps, 2, n2cFreqs)
    codonLLmax <- apply(codoncounts, 2, siteLLmax, pimatrix=codonProps)
    nucLLmax <- apply(nucleotidecounts, 2, siteLLmax, pimatrix=nucProps)
    meaNucID <- vapply(seq(1,3), function(i) seq(i,length(nucLLmax),by=3),
        FUN.VALUE=vector("numeric",length(nucLLmax)/3))
    nucLLsum <- apply(meaNucID, 1, function(i) sum(nucLLmax[i]))
    codonLLprob <- softmax(codonLLmax, TRUE)
    nucLLprob <- softmax(nucLLsum, TRUE)
    siteBF <- codonLLprob / nucLLprob
    red <- which( sum( colSums(codoncounts) > 0) == 1)
    a1out <- new("saturateBF", bf.pv=siteBF, nuLL=codonLLmax,
        altLL=nucLLsum, noInfo=red, nsites=ncol(codoncounts))
    return( a1out )
}

scanapp <- function(approach){
    meth <- toupper(approach)
    wrd <- "\nIncorrect `approach` specified. Expects \U2208 {`A`,`B`,`C`}"
    if(!(meth %in% c("A","B","C"))) stop(wrd, call.=FALSE)
    return(meth)
}

seqSaturation <- function(fastafile, approach, wsize, wless1, weightmax){
    method <- scanapp(approach)
    sanalyses <- switch(method,
        A = seqsat1(fastafile, wsize, wless1, weightmax),
        B = stop("\nApproach B is yet to be finalised."),
        C = stop("\nApproach C is yet to be finalised."))
    return(sanalyses)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

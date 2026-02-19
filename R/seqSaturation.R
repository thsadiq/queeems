# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                        Saturation Indices                        ><>< #
# ><>< ================================================================ ><>< #

betaspec <- function(nclass, less1prop, maxvalue){
    wrd1 <- "\nNumber of \U03B2 categories must be strictly > 3."
    wrd2 <- paste0("\nInappropriate \U03B2 < 1 proportion input:\n",
        "Must be such that at least two \U03B2 \U2264 1.")
    wrd3 <- paste0("\nIncorrect \U03B2 upper bound input:",
                "\nRequire 'maxvalue' > 1.")
    if(nclass <= 3) stop(wrd1, call.=FALSE)
    if(maxvalue <= 1) stop(wrd3, call.=FALSE)
    if( floor(less1prop * nclass) < 2) stop(wrd2, call.=FALSE)
    return( 0 )
}

fubarsplit <- function(nclass, less1prop, maxvalue){
    preamble <- betaspec(nclass, less1prop, maxvalue)
    fkappa <- floor(less1prop * nclass)
    fdeltatop <- (maxvalue - 1)^(1/3)
    fdeltabot <- nclass - fkappa
    fdelta <- fdeltatop / fdeltabot
    class1 <- seq(1, fkappa) / fkappa
    class2 <- 1 + ((seq(1, fdeltabot) * fdelta)^3)
    betavector <- c(class1, class2)
    return(betavector)
}

# oneBetaLL <- function(gdistBeta, baseCountMatrix){
#     basesize <- nrow(baseCountMatrix)
#     dirichletAlpha <- rgamma(basesize, 1, gdistBeta)
#     multinomialPI <- gtools::rdirichlet(1, dirichletAlpha)
#     siteLL <- apply(baseCountMatrix, 2, mnomLogl, pivector=multinomialPI)
#     return(siteLL)
# }

# betaSpaceLL <- function(Gnclass, Gless1, Gmaxv, baseCountMatrix){
#     allBetas <- fubarsplit(Gnclass, Gless1, Gmaxv)
#     sitellfull <- sapply(allBetas, oneBetaLL, baseCountMatrix=baseCountMatrix)
#     procSiteLL <- apply(sitellfull, 1, max)
#     return( procSiteLL)
# }

codon2int <- function(codon){
    nucIDs <- c(A=1, C=2, G=3, T=4)
    codonxters <- unlist( strsplit(codon, ""))
    nuConvert <- nucIDs[codonxters]
    names(nuConvert) <- NULL
    nuConvert <- nuConvert[!is.na(nuConvert)]
    return(nuConvert)
}

codonFracs <- function(nucFrac, allcodonID){
    codOneWeights <- apply(allcodonID, 2, function(j) prod(nucFrac[j]))
    codOneFrac <- softmax(codOneWeights, TRUE)
    return(codOneFrac)
}

piProcessor <- function(Gnclass, Gless1, Gmaxv){
    allweights <- fubarsplit(Gnclass, Gless1, Gmaxv)
    nucIndices <- gtools::permutations(Gnclass, 4)
    nucWeights <- apply(nucIndices, 1, function(i) allweights[i])
    nucFreqs <- apply(nucWeights, 2, softmax, maxadj=TRUE)
    allcodonID <- sapply(senseCodon, codon2int)
    codonFreqs <- apply(nucFreqs, 2, codonFracs, allcodonID=allcodonID)
    answer <- list(nucpi=nucFreqs, codonpi=codonFreqs)
    return(answer)
}

oneSiteLL <- function(pivector, baseCountMatrix){ ###
    siteLL <- apply(baseCountMatrix, 2, mnomLogl, pivector=pivector)
    return(siteLL)
}

seqsat1 <- function(fastafile, Gnclass, Gless1, Gmaxv){
    codoncounts <- basecensus( baseFrequency(fastafile, "codon"))
    nucleotidecounts <- basecensus( baseFrequency(fastafile, "dna"))
    generateFreqs <- piProcessor(Gnclass, Gless1, Gmaxv)
    # codonLL <- piComboLL(Gnclass, Gless1, Gmaxv, codoncounts)
    # dna0ll <- piComboLL(Gnclass, Gless1, Gmaxv, nucleotidecounts)
    # codonLL <- betaSpaceLL(Gnclass, Gless1, Gmaxv, codoncounts)
    # dna0ll <- betaSpaceLL(Gnclass, Gless1, Gmaxv, nucleotidecounts)
    dnaID <- sapply(seq(1,3), function(i) seq(i,length(dna0ll),by=3))
    dnaLLnew <- apply(dnaID, 1, function(i) sum(dna0ll[i]))
    siteBF <- codonLL / dnaLLnew
    red <- which( sum( colSums(codoncounts) > 0) == 1)
    a1out <- new("saturateBF", bf.pv=siteBF, nuLL=dnaLLnew,
        altLL=codonLL, noInfo=red, nsites=ncol(codoncounts))
    return( a1out )
}

# seqSaturation <- function(fastafile, basename, stats="bayes", cmismatch=FALSE){
#     omit3 <- ifelse(basename=="dna", paste(omit3rd), "non-applicable")
#     gb1 <- head( seq(gblow, 1, length=gbhalf+1), gbhalf)
#     gb2 <- tail( seq(1, gbupp, length=gbhalf+1), gbhalf)
#     countedbases <- baseFrequency(fastafile, basename)
#     noclassFreqs <- basecensus( countedbases )
#     if(omit3rd) noclassFreqs <- omitspot3(noclassFreqs)
#     invariants <-  which( colSums(noclassFreqs>0) == 1)
#     nleaves <- nseqs( countedbases )
#     sitesize <- ncol(noclassFreqs)
#     nbases <- nrow(noclassFreqs)
#     h1vect <- apply(noclassFreqs, 2, dnaEntropy,
#         nleaf=NA, fml="logl", llparam=rep(1/nbases,nbases))
#     h0matrix <- sapply(c(gb1,1,gb2), pivecLoop, baseCountMatrix=noclassFreqs)
#     siteIndices <- rbind(pnull=rowSums(h0matrix), palt=h1vect)
#     indExit <- new("saturated", llvalues=siteIndices, noInfo=invariants,
#         delete3rd=omit3, protein=basename, nsites=sitesize)
#     return(indExit)
# }

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

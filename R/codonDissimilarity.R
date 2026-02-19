# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><    Codon Census Matrix Matched by Nucleotide (Dis)Similarity.    ><>< #
# ><>< ================================================================ ><>< #

noncensus <- function(sitebaseID, siteData, triLogic, senses){
    cohort <- senses[ triLogic[sitebaseID,] ]
    baseNcensus <- any(cohort %in% siteData)
    return(baseNcensus)
}

siteNcensus <- function(siteData, triLogic, senses){
    dissFreqs <- matrix(0, nrow=61, ncol=1, dimnames=list(senses))
    availbases <- which(senses %in% siteData)
    ncvector <- sapply(availbases, noncensus, siteData=siteData,
                    triLogic=triLogic, senses=senses)
    dissFreqs[availbases] <- as.numeric(ncvector)
    return(dissFreqs)
}

codonDissimilarity.BStringSet <- function(bstringset, diffUnit){
    triLogic <- diffNucBinary(diffUnit)
    codonsequences <- bstringCodons(bstringset)
    codoncounts <- apply(codonsequences, 2,
        siteNcensus, triLogic=triLogic, senses=senseCodon)
    rownames(codoncounts) <- senseCodon
    return(codoncounts)
}

setMethod("codonDissimilarity",
    signature(bstringset="BStringSet", diffUnit="numeric"),
    codonDissimilarity.BStringSet)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

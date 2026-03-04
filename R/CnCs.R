# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  (Non-)Synonymous Count Matrix.                  ><>< #
# ><>< ================================================================ ><>< #

nscodoncount <- function(wcodon, sitecodons, nsmatrix){
    codonID <- which(senseCodon %in% wcodon)
    nslogical <- nsmatrix[codonID,]
    truecodons <- senseCodon[nslogical]
    nscensus <- sum( sitecodons %in% truecodons )
    return(nscensus)
}

getinvary <- function(sitecodons){
    sitesenses <- sitecodons[sitecodons %in% senseCodon]
    invariant <- length(unique(sitesenses)) == 1
    return(invariant)
}

CnCs <- function(bstringset, synonym, diffceil){
    nonsynarray <- nonSynonymous(synonym, diffceil)
    wildstring <- Biostrings::consensusString(bstringset)
    wildseqs <- Biostrings::BStringSet(wildstring)
    codonprotein <- bstringCodons(bstringset)
    seqinvaries <- apply(codonprotein, 2, getinvary)
    invarysites <- which( c(seqinvaries))
    wildcodons <- bstringCodons(wildseqs)
    answer <- vapply(seq(1,length(wildcodons)), function(i){
        nscodoncount(wildcodons[i], codonprotein[,i], nonsynarray) },
        FUN.VALUE=vector("integer",1))
    names(answer) <- NULL
    nsname <- ifelse(synonym, "Synonymous", "Non-Synonymous")
    reqcounts <- new("cncsframe", nscensus=answer,
        wildseq=wildcodons, nstype=nsname,
        maxdiff=diffceil, invariants=invarysites)
    return(reqcounts)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

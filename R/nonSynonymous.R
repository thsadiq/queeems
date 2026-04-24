# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                 (Non-)Synonymous Logical Matrix.                 ><>< #
# ><>< ================================================================ ><>< #

maxdiff <- function(diffceil){
    one <- diffNucBinary(1)
    two <- diffNucBinary(2)
    zero <- diffNucBinary(0)
    three <- diffNucBinary(3)
    if(diffceil == 0) output <- zero
    if(diffceil == 1) output <- zero | one
    if(diffceil == 2) output <- zero | one | two
    if(diffceil == 3) output <- zero | one | two | three
    return(output)
}

nsynvect <- function(codon, synonym){
    orderedTriplets <- aminoacid[senseCodon]
    syncodons <- orderedTriplets %in% aminoacid[codon]
    if(synonym){ output <- syncodons }else{ output <- !syncodons }
    codonID <- which( names( orderedTriplets) == codon)
    output[codonID] <- FALSE
    return(output)
}

nonSynonymous <- function(synonym=FALSE, diffceil=1){
    errMsg <- "Invalid `diffceil` input. It can only be set as 0, 1, 2 or 3."
    nsynArray <- vapply(senseCodon, nsynvect,
        synonym=synonym, FUN.VALUE=vector("logical",61))
    if(!(diffceil %in% seq(0,3))) stop( errMsg )
    mainMatrix <- maxdiff(diffceil) & nsynArray
    return(mainMatrix)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

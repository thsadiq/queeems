# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><             Nucleotide Discrepancy Matrix Generator.             ><>< #
# ><>< ================================================================ ><>< #

diffNucleotide <- function(codon){
    differ <- vapply(senseCodon,
        function(x) unlist(strsplit(x,"")) !=  unlist(strsplit(codon,"")),
        FUN.VALUE=vector("logical",3))
    return( colSums(differ) )
}

diffNucBinary <- function(diffUnit){
    errMsg <- "Invalid `diffUnit` input. It can only be set as 0, 1, 2 or 3."
    if(!(diffUnit %in% seq(0,3))) stop( errMsg )
    diffLogic <- vapply(senseCodon, function(x) diffNucleotide(x) == diffUnit,
        FUN.VALUE=vector("logical",61))
    return(diffLogic)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

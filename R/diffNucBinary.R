# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><             Nucleotide Discrepancy Matrix Generator.             ><>< #
# ><>< ================================================================ ><>< #

diffNucleotide <- function(codon){
    differ <- sapply(senseCodon,
        function(x) unlist(strsplit(x,"")) !=  unlist(strsplit(codon,"")) )
    return( colSums(differ) )
}

diffNucBinary <- function(diffUnit){
    errMsg <- "Invalid `diffUnit` input. It can only be set as 0, 1, 2 or 3."
    if(!(diffUnit %in% seq(0,3))) stop( errMsg )
    diffLogic <- sapply(senseCodon, function(x) diffNucleotide(x) == diffUnit)
    return(diffLogic)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                Approximate Multinomial Likelihood                ><>< #
# ><>< ================================================================ ><>< #

accessinputs <- function(countvector, pivector){
    scream2 <- "\n'pivector' and ''countvector must be non-negative values."
    scream0 <- "\nLengths of 'countvector' and 'pivector' must be equal."
    scream1 <- "'pivector' standardised because \U03A3\U1D70B \U2260 1."
    if(length(countvector) != length(pivector)) stop(scream0, call.=FALSE)
    if(any(countvector<0) | any(pivector<0)) stop(scream2, call.=FALSE)
    if( abs(sum(pivector) - 1) > 1e-8) warning(scream1, call.=FALSE)
    newpi <- pivector / sum( pivector )
    return(newpi)
}

mnomLogl <- function(countvector, pivector){
    stdpi <- accessinputs(countvector, pivector)
    llunits <- countvector * log(stdpi)
    llclean <- llunits[stdpi > 1e-12]
    llvalue <- sum(llclean)
    return(llvalue)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

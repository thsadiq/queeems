# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         Entropy Indices.                         ><>< #
# ><>< ================================================================ ><>< #

validBases <- function(countvector, nzero){
    msg0 <- "\nEntries of `countvector` should be positive integers."
    naBases <- is.na(countvector)|is.nan(countvector)|is.infinite(countvector)
    if( any(naBases) | any(countvector<0)) stop(msg0, call. = FALSE)
    if(nzero) countvector <- countvector[countvector != 0]
    return(countvector)
}

siFreq <- function(countvector, nleaf, nzero){
    msg0 <- "\nTotal base count can not be greater than leaf count."
    if(nleaf < sum(countvector)) stop(msg0, call. = FALSE)
    obsBases <- validBases(countvector, nzero)
    baseFreq <- obsBases / nleaf
    return(baseFreq)
}

si.shannon <- function(countvector, nleaf){
    baseProps <- siFreq(countvector, nleaf, TRUE)
    shID <-  (-1) * sum(baseProps * log(baseProps))
    return(shID)
}

si.renyi <- function(countvector, nleaf, ralpha){
    msg0 <- "\nSpecified Renyi \U03B1 \U2209 (0,1)."
    if(is.na(ralpha) | is.null(ralpha))  stop(msg0, call. = FALSE)
    if((ralpha<=0) | (ralpha>=1)) stop(msg0, call. = FALSE)
    baseProps <- siFreq(countvector, nleaf, TRUE)
    rID <-  (1 / (1 - ralpha) ) * log( sum(baseProps^ralpha) )
    return(rID)
}

checkFml <- function(fml){
    newTech <- tolower(fml)
    invalid <- !(newTech %in% c("shannon", "renyi"))
    msg <- paste("Incorrect `fml` specification.\n\t",
        " Input one of `shannon` or `renyi`.", sep="")
    if(invalid) stop(msg, call.=FALSE)
    return(newTech)
}

entropyindex <- function(countvector, nleaf, fml="shannon", ralpha=0.1){
    approach <- checkFml(fml)
    divIndex <- switch(approach,
        shannon = si.shannon(countvector, nleaf),
        renyi = si.renyi(countvector, nleaf, ralpha))
    indexest <- ifelse(all(countvector==0), 0, divIndex)
    return(indexest)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><   Convert Biostrings::readBStringSet Object to Codon Sequences   ><>< #
# ><>< ================================================================ ><>< #

miniGEN <- function(nucInit, bases){
    glueDNA <- paste(bases[seq(nucInit,nucInit+2)], collapse="")
    return(glueDNA)
}

codonExtant <- function(dnaVector){
    terminal <- length(dnaVector)
    validLength <- (length(dnaVector) %% 3) == 0
    if(!validLength) warning("DNA seq. length not codon-fit.", call.=FALSE)
    terminal <- ifelse(validLength, terminal, (length(dnaVector) %/% 3) * 3)
    newBases <- vapply(seq(1,terminal,3),
        miniGEN, bases=dnaVector, FUN.VALUE=vector("character",1))
    return(newBases)
}

bstringCodons.BStringSet <- function(bstringset){
    newcseq <- t( apply( as.matrix(bstringset), 1, codonExtant))
    return(newcseq)
}

setMethod("bstringCodons", signature("BStringSet"), bstringCodons.BStringSet)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

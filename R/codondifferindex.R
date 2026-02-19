# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><      Molecular Entropy Generated with Codon Mismatch Counts      ><>< #
# ><>< ================================================================ ><>< #

codondifferindex <- function(fastafile, diffUnit, fml, ralpha=NA){
    if(!(diffUnit%in%seq(1,3))) stop("\nInput `diffUnit` \U2284 {1,2,3}")
    basestrings <- Biostrings::readBStringSet(fastafile)
    differoutput <- codonDissimilarity(basestrings, diffUnit)
    nouse <-  which(colSums( differoutput) == 0)
    nleaves <- length(basestrings)
    eindex <- apply(differoutput, 2,
        entropyindex, nleaf=nleaves, fml=fml, ralpha=ralpha)
    entropyexit <- new("siteindices",
        eindices=eindex, noinfo=nouse, meth=fml, alphaR=as.numeric(ralpha))
    return(entropyexit)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

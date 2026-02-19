# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><              Global (Non-)Synonymous Entropy Index.              ><>< #
# ><>< ================================================================ ><>< #

cncsentropy <- function(fastafile, synonym, diffceil, fml, ralpha=NA){
    basestrings <- Biostrings::readBStringSet(fastafile)
    cnsoutput <- CnCs(basestrings, synonym, diffceil)
    nleaves <- length(basestrings)
    entcount <- nsfreqs(cnsoutput)
    countsites <- length(entcount)
    basesize <- nleaves * countsites
    nosign <- entcount[entcount > 0]
    entvalue <- entropyindex(entcount, nleaf=basesize, fml=fml, ralpha=ralpha)
    centout <- new("geneindex", eindex=entvalue, noinfo=nosign,
        meth=fml, alphaR=as.numeric(ralpha),
        nbases=countsites, synonym=synonym)
    return(centout)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

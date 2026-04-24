# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><       Molecular Entropy Generated For Individual Base Site       ><>< #
# ><>< ================================================================ ><>< #

molentropy <- function(fastafile, basename="dna", fml="shannon", ralpha=NA){
    molesumary <- baseFrequency(fastafile, basename)
    molecounts <- basecensus( molesumary )
    leafprevs <- nseqs( molesumary )
    unibase <-  which( colSums(molecounts) == 1)
    tropyI <- apply(molecounts, 2, entropyindex,
        nleaf=leafprevs, fml=fml, ralpha=ralpha)
    exitbox <- new("siteindices", eindices=tropyI,
        noinfo=unibase, meth=fml, alphaR=as.numeric(ralpha))
    return(exitbox)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

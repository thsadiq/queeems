# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                 Rid Sequence of Saturated Sites.                 ><>< #
# ><>< ================================================================ ><>< #

siteIDcheck <- function(deadsites, nsites){
    rule1 <- any(deadsites <= 0)
    rule2 <- is.na( sum(deadsites))
    rule3 <- length(deadsites) == 0
    rule4 <- any(!(deadsites < Inf))
    text1 <- "\nElements of `deadsites` must be positive integers."
    text2 <- "\nShould be removing less than the total site size."
    if(rule1 | rule2 | rule3 | rule4) stop(text1, call.=FALSE)
    if(length(deadsites) > nsites) stop(text2, call.=FALSE)
    text3 <- "\nAt least one of `deadsites` is undefined."
    if(max(deadsites) > nsites) stop(text3, call.=FALSE)
    return(deadsites)
}

seqfilter <- function(fastafile, deadsites, basename, write=FALSE){
    cleanbase <- baseCheck(basename)
    baseseqs <- switch(cleanbase,
        aa = as.matrix( Biostrings::readAAStringSet(fastafile)),
        dna = as.matrix( Biostrings::readDNAStringSet(fastafile)),
        codon = bstringCodons( Biostrings::readBStringSet(fastafile)))
    cleansiteid <- siteIDcheck(deadsites, ncol(baseseqs))
    cleanmat <- baseseqs[,-cleansiteid]
    genetext <- apply(cleanmat, 1, paste, collapse="")
    editedseqs <- switch(cleanbase,
        aa = Biostrings::AAStringSet(genetext),
        dna = Biostrings::DNAStringSet(genetext),
        codon = Biostrings::BStringSet(genetext))
    if(write != FALSE){
        filepath <- as.character(write)
        Biostrings::writeXStringSet(editedseqs, filepath)
        message("\n\nUpdated sequence file saved as:\n\t", filepath, "\n")
    }
    return(editedseqs)
}




# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

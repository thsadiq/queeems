# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  Methods for `cncsframe` Class.                  ><>< #
# ><>< ================================================================ ><>< #

setMethod("show", "cncsframe", function(object) {
    cat("\n# :::::::")
    cat("\nCounts generator:\t\t ", object@nstype, "codons")
    cat("\nNon-informative codon sites:\t ", sum(object@nscensus == 0),
        paste0("(of ", length(object@wildseq), ")"))
    cat("\nMax. nuc. mismatch tolerance:\t ", object@maxdiff)
    cat("\n# :::::::\n\n")
    }
)
setMethod("nORs", "cncsframe", function(cNS) cNS@nstype )
setMethod("nsfreqs", "cncsframe", function(cNS) cNS@nscensus )
setMethod("nucbalance", "cncsframe", function(cNS) cNS@maxdiff )
setMethod("wildtriplets", "cncsframe", function(cNS) cNS@wildseq )
setMethod("nonvaries", "cncsframe", function(cNS) cNS@invariants )

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

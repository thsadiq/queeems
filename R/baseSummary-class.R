# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                 Methods for `baseSummary` Class.                 ><>< #
# ><>< ================================================================ ><>< #

setMethod("show", "baseSummary", function(object) {
    cat("\n# :::::::")
    temp <- as.numeric( table( colSums(object@seqsNumerics)))
    cat("\nProtein base type:\t\t ", object@baseType)
    cat("\nNumber of sequences:\t\t ", object@extantSize)
    cat(paste("\nDimension of numeric output:\t", sep=""))
    cat(paste("(base_size=", nrow(object@seqsNumerics), ")", sep=""))
    cat(paste("-by-(site_size=", ncol(object@seqsNumerics), ")\n", sep=""))
    cat("# :::::::\n\n")
    }
)
setMethod("nseqs", "baseSummary", function(bsumm) bsumm@extantSize )
setMethod("datatype", "baseSummary", function(bsumm) bsumm@baseType )
setMethod("basecensus", "baseSummary", function(bsumm) bsumm@seqsNumerics )

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

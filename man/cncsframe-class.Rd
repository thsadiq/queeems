\name{cncsframe-class}
\docType{class}
\alias{cncsframe}
\alias{cncsframe-class}
\alias{show, cncsframe-method}
\alias{nORs, cncsframe-method}
\alias{nsfreqs, cncsframe-method}
\alias{nonvaries, cncsframe-method}
\alias{nucbalance, cncsframe-method}
\alias{wildtriplets, cncsframe-method}

\title{Synonymous and Non-Synonymous Frequencies}

\description{
    Facilitate smooth extraction and utility of outputs from synonymous
    and non-synonymous count analyses of molecular protein data.
}

\section{Objects from the Class}{
    Objects of this (\code{cncsframe}) class can be created by calls of
    the form \code{new("cncsframe", nscensus, wildseq, nstype, maxdiff,
    invariants)}.
}

\section{Slots}{
    \describe{
        \item{\code{nscensus}:}{site-wise counts of codons that are synonymous
            or non-synonymous to the wild-type codon.}
        \item{\code{wildseq}:}{sequence that comprise of the most frequent
            codon at each site.}
        \item{\code{nstype}:}{indication of the type, i.e. \code{synonymous}
            or \code{non-synonymous}, of counts returned.}
        \item{\code{maxdiff}:}{maximum number of nucleotide mismatch
            that is tolerated.}
        \item{\code{invariants}:}{indices of the non-informative codon sites.}
    }
}

\section{Methods}{
    \describe{
        \item{show:}{\code{signature(object="cncsframe")} snapshot of
            output features.}
        \item{nORs:}{\code{signature(cNS="cncsframe")} computational
            approach, `\code{synonymous}` or `\code{non-synonymous}`.}
        \item{nsfreqs:}{\code{signature(cNS="cncsframe")} site-wise
            counts of codons that are (non-)synonymous to the wild-type.}
        \item{nucbalance:}{\code{signature(cNS="cncsframe")} maximum
            number of nucleotide mismatch tolerated per codon pair.}
        \item{wildtriplets:}{\code{signature(cNS="cncsframe")} a vector
            of the identified wild-type codon for each  site of the input
            protein data.}
        \item{nonvaries:}{\code{signature(cNS="cncsframe")} locations
            of sites that were detected to be non-informative.}
    }
}

\section{Details}{
    Consult the \code{\link{CnCs}} function for more elaborate descriptions
    of the slots and method outputs.
}

\author{
    Hassan Sadiq
}

\seealso{
    The \code{\link{CnCs}} function that returns \code{cncsframe} objects.
}

\examples{
    testsq <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
                "CCTCAGATCACTCTTGTCTCAATAAAAGTA",
                "TGG-AGCGACCCCTTGTCTCAATAAAAATA")
    codonarray <- Biostrings::BStringSet(testsq)
    tested <- CnCs(codonarray, FALSE, 3)
    wildtriplets(tested)
    nsfreqs(tested)
    tested
}

\keyword{classes}

\name{CnCs}
\alias{CnCs}

\title{Synonymous and Non-Synonymous Codon Match Counts}

\description{
    Crude counts of the number of synonymous or non-synonymous matches at
    each codon site of input empirical molecular data.
}

\usage{
    CnCs(bstringset, synonym, diffceil)
}

\arguments{
    \item{bstringset}{Biostrings::Bstringset object that carries an
        empirical protein sequence.}
    \item{synonym}{logical entry that indicates whether synonymous
        [\code{TRUE}] or non-synonymous [\code{FALSE}] counts are needed.}
    \item{diffceil}{maximum number of nucleotide mismatches that
        is tolerated.}
}

\details{
    Each element of the \code{nscensus} output provides information about the
    respective codon site. Consider the case where \code{synonym = TRUE}. At
    each site, the most frequent codon is identified as the wildtype base.
    The number of spots at the site that are filled by separate but synonymous
    codons to the wild-type is recorded. Ambiguities and/or missing
    observations are ignored. Sites with zero counts are considered
    non-informative. Observe that such sites are not necessarily invariant.
    For example, a site with three unique non-synonymous codon bases will
    only be informative when non-synonymous counts are sought, i.e. if
    \code{synonym = FALSE}. The \code{diffceil} input implies the maximum,
    not the exact, number of nucleotide mismatches allowed. For example,
    \code{diffceil = 2} accommodates zero, one or two nucleotide mismatches.
}

\value{
    A \code{\link{cncsframe}} object that contains the following.
    \item{nscensus}{an array of integers that may be accessed with the
        \code{\link{nsfreqs}} function. They represent site-wise counts,
        compared to the identified wild-type, of [non-]synonymous
        occurrences.}
    \item{wildseq}{an array of the most frequent, i.e. the identified
        wild-type, codon at each site of the input molecular data.}
    \item{nstype}{a string that describes whether the entries of
        \code{nscensus} were obtained with respect to synonymous or
        non-synonymous counts.}
    \item{maxdiff}{maximum number of nucleotide mismatches tolerated when
        comparing to the wild-type.}
    \item{invariants}{indices of invariant codon sites, after ambiguities
        and missingness have been accommodated.}
}

\author{
    Hassan Sadiq
}

\seealso{
    The \code{\link{cncsframe}} class that contains additional informationn
    about how to access the output object. The \code{\link{nonSynonymous}}
    and \code{\link{codonDissimilarity}} functions may be useful resources
    for clarifications about the \code{maxdiff} entry.
}

\examples{
    indata <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
                "CCTCAGATCACTCTTGTCTCAATAAAAGTA",
                "TGG-AGCGACCCCTTGTCTCAATAAAAATA")
    codonarray <- Biostrings::BStringSet(indata)
    tested <- CnCs(codonarray, FALSE, 3)
    tested
}

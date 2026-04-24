\name{cncsentropy}
\alias{cncsentropy}

\title{Global [Non-]Synonymous Entropy Index for Protein Sequence}

\description{
    Obtain homogeneous entropy estimate for a set of homologous protein
    sequences, based on observed synonymous or non-synonymous codon
    diversities.
}

\usage{
    cncsentropy(fastafile, synonym=FALSE, diffceil=1, fml="shannon", ralpha=NA)
}

\arguments{
    \item{fastafile}{path to a FASTA format file that contains a protein
        sequence.}
    \item{synonym}{logical entry that indicates whether synonymous
        [\code{TRUE}] or non-synonymous [\code{FALSE} (default)] counts
        are preferred.}
    \item{diffceil}{maximum number of nucleotide mismatches that
        is tolerated. Default value is \code{1}.}
    \item{fml}{character input that can be one of \code{"shannon"} (default)
        or \code{"renyi"}.}
    \item{ralpha}{a fraction between 0 and 1 that is needed for Renyi
        index calculation.}
}

\details{
    Entropy calculations are largely functions of relative counts. See the
    \code{\link{entropyindex}} function for the mathematical expressions
    and some useful citations. Obtaining site-specific entropy values,
    where the proportions represent the relative occurrence of the
    protein bases per site, is not ambiguous. With respect to obtaining
    a single entropy estimate for an entire alignment in terms of
    synonymous and non-synonymous relationships, the process is less
    straightforward. The proportions used for the calculations were
    obtained as follows.
    \deqn{
        p^{\code{syn}}_{i} = \frac{ n^{\code{syn}}_{i} }
        {\sum_{i=1}^{k} n^{\code{syn}}_{i}},
    }
    where \eqn{p^{\code{syn}}_{i}} represents the proportion of observed
    synonymous substitutions at site \eqn{i = 1, 2, \ldots, k} and
    \eqn{n^{\code{syn}}_{i}} is the number of observed synonymous
    substitutions with reference to the most frequently occurring codon
    at the site. For example, the expression for the corresponding
    synonymous Shannon entropy for a protein sequence alignment then
    becomes,
    \deqn{
        H^{\code{syn}}_{} = -\sum_{i=1}^{k} p^{\code{syn}}_{i} \times
            \log p^{\code{syn}}_{i}.
    }
    Similar explanation holds for estimates obtained with respect to
    non-synonymous substitutions or from using the Renyi approach.
}

\value{
    A \code{\link{geneindex}} object.
}

\author{
    Hassan Sadiq
}

\seealso{
    For further clarification about the operation of the \code{cncsentropy}
    function and extensive understanding of the interpretation of its output,
    consult the constituent \code{\link{CnCs}} and \code{\link{entropyindex}}
    functions. There are also the \code{\link{codondifferindex}} and
    \code{\link{molentropy}} functions that are capable of generating
    site-wise indices without distinguishing between synonymous or
    non-synonymous codon relationships.
}

\references{
    Shannon, C. E. (1948). A Mathematical Theory of Communication, \emph{The
    Bell System Technical Journal} \bold{27}(3): 379-423.
}

\examples{
    sqdatm <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
                "CCTCAGATCACTCTTGTCTCAATAAAAGTA",
                "CCC-AGCGACCCCTTGTCTCAATAAAAATA")
    codonData <- Biostrings::BStringSet(sqdatm)
    names(codonData) <- c("s1","s2","s3")
    testdna <- file.path(tempdir(), "testDna.fasta")
    Biostrings::writeXStringSet(codonData, testdna)
    testing <- cncsentropy(testdna, TRUE, 2, "shannon")
    print(testing)
}

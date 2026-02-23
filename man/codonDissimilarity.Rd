\name{codonDissimilarity}
\alias{codonDissimilarity}
\alias{codonDissimilarity.BStringSet}

\title{Codon Census Matrix Matched by Nucleotide Dissimilarity}

\description{
    Generates a matrix that contains numerical details of the distribution
    of bases at each codon position in a protein sequence.
}

\usage{
    codonDissimilarity(bstringset, diffUnit)
}

\arguments{
    \item{bstringset}{A BStringSet object.}
    \item{diffUnit}{An integer between the 0 and 3 interval.}
}

\details{
    The input \code{BStringSet} object can be created with the
    \code{\link{Biostrings::readBStringSet}} function. \code{diffUnit}
    represents the number of nucleotide mismatches that is preferred. Let
    the unique codons at each genetic sequence alignment site be denoted
    as \code{i = 1, 2, ..., k}, where \code{k} can not be more than the
    number of sequences that make up the data. For each \code{i} at the
    different sites, a value of 1 is recored if any of the other codons
    at the site differ from codon \code{i} by exactly \code{diffUnit}
    nucleotide positions. The process is repeated independently for each
    site and the outcome is entered into the corresponding column of the
    output matrix. Note the special case of \bold{\code{diffUnit=0}}. The
    output will be a binary matrix. Row \code{r} column \code{c} of the
    matrix will contain a 0 if codon \code{r} was observed at site
    \code{c} of the input protein sequence.
}

\value{
    A 61-by-\code{nsites} dissimilarity matrix that captures the site-wise
    distribution of codon bases evident from the input protein data, where
    \code{nsites} is the number of codon sites. As explained as part of
    \bold{Details} above, codon dissimilarity is determined by specified
    degree of nucleotide mismatch.
}

\seealso{
    \code{\link{queeemsExtdata}} that returns paths to the example
    external DNA data that can be processed with the described function.
    \code{Biostrings::readBStringSet}, a necessary feeder function may
    also be a useful resource. The \code{\link{baseFrequency}} function
    implicitly uses this function.
}

\author{
    Hassan Sadiq
}

\examples{
    dataPath <- queeemsExtdata("PI.fasta")
    bstringdata <- readBStringSet(dataPath)
    testing <- codonDissimilarity(bstringdata, 1)
    head(testing)[,seq(1,6)]
}

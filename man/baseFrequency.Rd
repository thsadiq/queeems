\name{baseFrequency}
\alias{baseFrequency}

\title{Counts Indicating Protein Base Prevalence}

\description{
    Obtain, from empirical molecular data, frequencies of each base (DNA,
    codon or amino acid) per data site.
}

\usage{
    baseFrequency(fastafile, basename="dna")
}

\arguments{
    \item{fastafile}{path to a FASTA format file that contains a protein
        sequence.}
    \item{basename}{identification of the protein constituent base.}
}

\details{
    \code{basename} should be a character input and one of `dna` (defualt),
    `codon` or `aa`. Any other input will trigger an error. A warning message
    may also be returned if there is an indication of sequence-basename
    mismatch. For example if `dna` is specified for a protein sequence
    saved in terms of amino acid bases.
}

\value{
    A \code{\link{baseSummary}} object that contains a matrix of site-wise
    base frequencies that are labelled with the corresponding base notation.
}

\author{
    Hassan Sadiq
}

\seealso{
    The \code{\link{baseSummary}} class for further description of the
    output and other similar fuctions such as \code{\link{CnCs}} and
    \code{\link{diffNucBinary}}. Sample FASTA files may also be accessed
    using the \code{\link{queeemsExtdata}} function.
}

\examples{
    testData <- queeemsExtdata("PI.fasta")
    tested <- baseFrequency(testData, "codon")
    tested
}

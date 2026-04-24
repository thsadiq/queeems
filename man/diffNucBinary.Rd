\name{diffNucBinary}
\alias{diffNucBinary}

\title{Compare Codon Pairs to Identify Nucleotide Mismatch}

\description{
    Generate a binary matrix that identifies pairs of sense codons that differ
    by a specified degree.
}

\usage{
    diffNucBinary(diffUnit=1)
}

\arguments{
    \item{diffUnit}{A non-negative integer. Defaults to 1.}
}

\details{
    \code{diffUnit} should ideally be a non-negative value that is less than
    4. The codons are matched by nucleotide bases and their positions. For
    example, comapring codon \code{`ACA`} and codon \code{`CAA`} will return
    \code{TRUE} only when \code{diffUnit = 2} in accordance with the
    discrepancies at the first and the second nuleotide positions. Observe
    that when \code{diffUnit = 0}, a logical diagonal matrix is returned.
}

\value{
    A 61-by-61 logical matrix. \code{TRUE} will appear at the row and column
    intersection of the sense codons that differ by the number of nucleotide
    bases specified by \code{diffUnit}. \code{FALSE} will appear otherwise.
}

\author{
    Hassan Sadiq
}

\examples{
    test1 <- diffNucBinary(1)
    head(test1)[,seq(1,6)]
}

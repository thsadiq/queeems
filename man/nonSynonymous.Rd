\name{nonSynonymous}
\alias{nonSynonymous}

\title{Identify Synonymous and Non-Synonymous Codons}

\description{
    Create a 61-by-61 symmetric logical matrix that indicates the sense
    codons that are synonymous or non-synonymous to the corresponding
    row or column codon.
}

\usage{
    nonSynonymous(synonym, diffceil)
}

\arguments{
    \item{synonym}{logical input that indicates whether synonymous [TRUE]
        or non-synonymous identifier is required.}
    \item{diffceil}{an interger that states the preferred maximum number
        of nucleotide mismatch.}
}

\details{
    \code{diffceil} should ideally be a non-negative value that is less than
    4. The codons are matched by nucleotide bases and their positions. For
    example, when \code{diffceil = 2} and \code{synonym = FALSE}, the output
    matrix column associated with codon `\code{AAA}` will contain \code{TRUE}
    for all codons that are non-synonymous to `\code{AAA}` and that differ
    from `\code{AAA}` at not more than two nucleotide positions. Codons
    that do not satisfy that criterion will be assigned \code{FALSE}. The
    same interpretation applies to the respective rows of the output matrix
    corresponding to the other sense codons. The output matrix is symmetric,
    so all the interpretations hold similarly row-wise. The rows and columns
    of the output matrix are arranged and labelled alphabetically with
    respect to the IUPAC nucleotide triplet naming style.
}

\value{
   A 61-by-61 logical matrix. See \bold{Details}.
}

\author{
    Hassan Sadiq
}

\seealso{
    Explanations that accompany the \code{\link{diffNucBinary}} function
    might provide better operational understanding.
}

\examples{
    testing <- nonSynonymous(FALSE, 2)
    head(testing)[,seq(1,6)]
}

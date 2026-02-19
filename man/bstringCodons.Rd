\name{bstringCodons}
\alias{bstringCodons}
\alias{bstringCodons.BStringSet}

\title{Convert a BStringSet Object to Codon Sequences}

\description{
    Converts a Biostrings::BStringSet object to a matrix of codons.
}

\usage{
    bstringCodons(bstringset)
}

\arguments{
    \item{bstringset}{A BStringSet object.}
}

\value{
    A matrix of codon bases. Each column of the matrix will contain the site
    data and each row will hold the respective sequence that make up the input
    protein data.
}

\author{
    Hassan Sadiq
}

\examples{
    dataPath <- queeemsExtdata("PI.fasta")
    bstringdata <- readBStringSet(dataPath)
    testing <- bstringCodons(bstringdata)
    head(testing)[,seq(1,6)]
}

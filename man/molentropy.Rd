\name{molentropy}
\alias{molentropy}

\title{Molecular Entropy Generated Independently per Base Site}

\description{
    Generate information entropy estimate independently for protein
    sequence sites.
}

\usage{
    molentropy(fastafile, basename, fml, ralpha=NA)
}

\arguments{
    \item{fastafile}{path to a FASTA format file that contains a protein
        sequence.}
    \item{basename}{specification of how the data must be encoded. It can
        be \code{'dna'}, \code{'codon'} or \code{'aa'} for DNA, codon or
        amino acid, respectively.}
    \item{fml}{character input that can be one of \code{"shannon"} or
        \code{"renyi"}.}
    \item{ralpha}{a fraction between 0 and 1 that is needed for Renyi
        index calculation.}
}

\value{
    A \code{\link{siteindices}} object.
}

\author{
    Hassan Sadiq
}

\seealso{
    For further clarification about the operation of the \code{molentropy}
    function and extensive understanding of the interpretation of its output,
    consult the \code{\link{baseFrequency}} and \code{\link{entropyindex}}
    functions. There is also the \code{\link{codondifferindex}} function
    that can be used to generate site-wise entropy indices and the
    \code{\link{cncsentropy}} function that is suitable for homogeneous
    entropy estimation.
}

\examples{
    sqdatm <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
                "CCTCAGATCACTCTTGTCTCAATAAAAGTA",
                "CCC-AGCGACCCCTTGTCTCAATAAAAATA")
    moleculeDta <- Biostrings::BStringSet(sqdatm)
    names(moleculeDta) <- c("s1","s2","s3")
    testprotein <- file.path(tempdir(), "testprotein.fasta")
    Biostrings::writeXStringSet(moleculeDta, testprotein)
    test1 <- molentropy(testprotein, "dna", "renyi", .35)
    test2 <- molentropy(testprotein, "codon", "shannon")
    print(test1)
    print(test2)
}

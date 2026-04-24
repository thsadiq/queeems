\name{codondifferindex}
\alias{codondifferindex}

\title{Molecular Entropy Estimated from Codon Mismatches}

\description{
    Obtain information entropy index for each codon site in a protein
    sequence, as a function of mismatches observed at nucleotide positions.
}

\usage{
    codondifferindex(fastafile, diffUnit=1, fml="shannon", ralpha=NA)
}

\arguments{
    \item{fastafile}{path to a FASTA format file that contains a protein
        sequence.}
    \item{diffUnit}{number of nucleotide mismatches between codon pairs
        to use for the calculations.}
    \item{fml}{character input that can be one of \code{"shannon"} (default)
        or \code{"renyi"}.}
    \item{ralpha}{a fraction between 0 and 1 that is needed for Renyi
        index calculation.}
}

\value{
    An object of class \code{\link{siteindices}}.
}

\author{
    Hassan Sadiq
}

\seealso{
    This function is a combination of the operations of the
    \code{\link{codonDissimilarity}} and the \code{\link{entropyindex}}
    functions. The two referenced functions contain further discussions
    about the \code{diffUnit}, \code{fml} and \code{ralpha} inputs. See
    the \code{\link{siteindices}} class for an extensive detail about
    how to work with the output object. The \code{\link{molentropy}}
    function is an example of an alternative approach of estimating
    site-wise entropy indices accessible within the \code{\link{queeems}}
    package. 
}

\examples{
    proteinFile <- queeemsExtdata("VertCOI.fasta")
    testEG <- codondifferindex(proteinFile, 1, "shannon")
    summary(testEG)
    print(testEG)
}

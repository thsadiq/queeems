\name{n2cFreqs}
\alias{n2cFreqs}

\title{Generate Codon Frequencies from Nucleotide Frequencies}

\description{
    Obtain frequencies of sense codons as a function of the frequencies of
    nucleotide bases
}

\usage{
    n2cFreqs(nucFracs)
}

\arguments{
    \item{nucFracs}{vector of nucleotide frequencies.}
}

\details{
    This provides a straightforward way to generate codon frequencies as
    described for the popular F3x4 model. For a given vector (\code{nucFracs}
    = \eqn{\boldsymbol{\pi}^{(n)}} = \eqn{\{\pi_{\code{A}}^{(n)},
    \pi_{\code{C}}^{(n)}, \pi_{\code{G}}^{(n)}, \pi_{\code{T}}^{(n)}\}}) of
    nucleotide frequencies, each non-standardised sense codon frequency,
    \eqn{\pi_{i}^{(c)}} for \code{i = 1, 2, ..., 61} is obtained as the product
    the product of the nucleotide triplet that form the codon. For example, for
    \eqn{\pi_{i}^{(c)}} = \eqn{\pi_{\code{AGT}}^{(c)}}, the unstandardised
    codon frequency is estimated as \eqn{ \pi_{\code{A}}^{(n)} \times
    \pi_{\code{G}}^{(n)} \times \pi_{\code{T}}^{(n)} }. Please note that, if
    the contents of the input vector are not named, the elements are assigned
    base names (\code{A}, \code{C}, \code{G}, \code{T}) in order of appearance
    in the input.
}

\value{
    A numeric vector of codon frequencies that are standardised to sum to one,
    ordered and named in terms of the IUPAC nucleotide triplet nomenclature.
}

\author{
    Hassan Sadiq
}

\references{
    Goldman, N. and Yang, Z. (1994) A Codon-Based Model of Nucleotide
    Substitution for Protein-Coding DNA Sequences, \emph{Mol. Biol.
    Evol.} \bold{11}(5): 725-736.
}

\examples{
    codonFrac1 <- n2cFreqs( c(A=0.3, C=0.1, G=0.4, T=0.2) )
    codonFrac2 <- n2cFreqs( c(0.3, 0.1, 0.4, 0.2) )
    print(codonFrac1)
    print(codonFrac2)
}

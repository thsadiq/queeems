\name{baseSummary-class}
\docType{class}
\alias{baseSummary}
\alias{baseSummary-class}
\alias{nseqs}
\alias{datatype}
\alias{basecensus}
\alias{show,baseSummary-method}
\alias{nseqs,baseSummary-method}
\alias{datatype,baseSummary-method}
\alias{basecensus,baseSummary-method}

\title{Protein Base Composition Numerics}

\description{
    Ensure that numerical extracts from protein data are stored and handled
    elegantly.
}

\section{Objects from the Class}{
    Objects of this (\code{baseSummary}) class can be created by calls of
    the form \code{new("baseSummary", baseType, extantSize, seqsNumerics)}.
}

\section{Slots}{
    \describe{
        \item{\code{baseType}:}{description of the bases that make up the
            corrsponding protein sequence.}
        \item{\code{extantSize}:}{number of sequences in the protein data.}
        \item{\code{seqsNumerics}:}{numerical data that describe the
            distribution of the primary bases in the captured protein
            sequences.}
    }
}

\section{Methods}{
    \describe{
        \item{show:}{\code{signature(object="baseSummary")} useful output
            highlights.}
        \item{nseqs:}{\code{signature(bsumm="baseSummary")} number of
            protein sequences.}
        \item{datatype:}{\code{signature(bsumm="baseSummary")} how protein
            data positions must be interpreted.}
        \item{basecensus:}{\code{signature(bsumm="baseSummary")} numerical
            summary of the empirical protein base distribution.}
    }
}

\section{Details}{
    The \code{baseType} slot should be one of \code{`aa`}, \code{`codon`}
    or \code{`dna`} depending on whether the protein sequence sites must
    be analysed as `amino acid`, `codon` or `DNA` positions. The rows of the
    element of the \code{seqsNumerics} slot are defined with respect to the
    alphabetically ordered IUPAC abbreviations of the corresponding
    \code{baseType} input. Each column of the matrix represents the relevant
    protein site.
}

\author{
    Hassan Sadiq
}

\seealso{
    The \code{\link{baseFrequency}} function that returns \code{baseSummary}
    objects.
}

\examples{
    proteinFile <- queeemsExtdata("PI.fasta")
    testEG <- baseFrequency(proteinFile, "codon")
    head( basecensus(testEG))[,seq(1,6)]
    datatype(testEG)
    nseqs(testEG)
    testEG
}

\keyword{classes}

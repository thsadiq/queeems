\name{seqfilter}
\alias{seqfilter}

\title{Delete Subset of Molecular Sequences}

\description{
    Extract targetted base sites from molecular sequence data.
}

\usage{
    seqfilter(fastafile, deadsites, basename, write=FALSE)
}

\arguments{
    \item{fastafile}{full path to the file that contains the sequence that
        needs to be subsetted.}
    \item{deadsites}{integer indices of the base locations that should be
        removed from the sequence in \code{fastafile}.}
    \item{basename}{one of \code{`dna`}, \code{`codon`} or \code{`aa`} that
        indicates how the sequence data should be captured.}
    \item{write}{optional full path to a file where the output sequence
        should be saved.}
}

\details{
    If entered, \code{write} should be the full path to the file where the
    data is expected to be saved. The input is not always interpreted as
    logical. If \code{write = FALSE} is entered (default), no file will be
    saved. However, \code{write = TRUE} will save the edited sequence output
    in a file named \code{`TRUE`} in the working directory.
}

\value{
    The extracted molecular sequence data as a \code{Biostrings::BStringSet}
    object. If \code{write} is specified, the edited sequence is also saved
    in the filepath stated.
}

\author{
    Hassan Sadiq
}

\examples{
    egdata <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
                "CCTCAGATCACTCTTGTCTCAATAAAAGTA",
                "TGG-AGCGACCCCTTGTCTCAATAAAAATA")
    dnaseqs <- Biostrings::DNAStringSet(egdata)
    names(dnaseqs) <- paste0("seq",seq(1,3))
    egpath <- file.path(tempdir(), "egdata.fas")
    Biostrings::writeXStringSet(dnaseqs, egpath)
    cleanseqs <- seqfilter(egpath, c(4,6,8), "dna")
    print(cleanseqs)
}

\name{saturatePV-class}
\docType{class}
\alias{saturatePV}
\alias{saturatePV-class}
\alias{PVs, saturatePV-method}
\alias{show, saturatePV-method}
\alias{summary, saturatePV-method}

\title{Substitution Saturation Indices: Frequentist Paradigm}

\description{
    Allow for neat and efficient handling of outputs with respect
    to substitution saturation indices calculations following the
    frequentist principles.
}

\section{Objects from the Class}{
    Creating \code{saturatePV} objects may be achieved with calls of the
    form: \code{new("saturatePV", bf.pv, nuLL, altLL, noInfo, nsites)}.
}

\section{Slots}{
    Accessible slots of \code{saturatePV} objects are the same as the
    slots of the \code{\link{saturateBF}} class object from where it
    inherits, note the updated description of the following.
    \describe{
        \item{\code{bf.pv}:}{p-values obtained from independent
            implementation of Chi-squared likelihood ratio test
            for each base site.}
    }
}

\section{Methods}{
    \describe{
        \item{show:}{\code{signature(object="saturatePV")} output highlights
            that are self-explanatory.}
        \item{summary:}{\code{signature(object="saturatePV")} interesting
            saturation statistics that are informed by Chi-squared likelihood
            ratio test.}
        \item{PVs:}{\code{signature(satupe="saturatePV")} p-values recorded
            for each base site.}
    }
}

\section{Details}{
    Unlike for objects of the inherited \code{\link{saturateBF}} class,
    p-values are returned as the measure of saturation significance for
    each base site.
}

\author{
    Hassan Sadiq
}

\seealso{
    The inherited \code{\link{saturateBF}} class and the
    \code{\link{seqSaturation}} function that produces objects of this
    class for more information.
}

\examples{
    moleFile <- queeemsExtdata("II.fasta")
    assess <- seqSaturation(moleFile, 3)
    assess
}

\keyword{classes}

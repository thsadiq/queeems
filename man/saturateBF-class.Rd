\name{saturateBF-class}
\docType{class}
\alias{saturateBF}
\alias{saturateBF-class}
\alias{BFs, saturateBF-method}
\alias{show, saturateBF-method}
\alias{LogL0, saturateBF-method}
\alias{LogL1, saturateBF-method}
\alias{summary, saturateBF-method}
\alias{sitecount, saturateBF-method}
\alias{nonvaries, saturateBF-method}

\title{Substitution Saturation Output: Bayesian Paradigm}

\description{
    Allow for neat and efficient handling of outputs from Bayesian
    site-wise substitution saturation analyses.
}

\section{Objects from the Class}{
    Creating \code{saturateBF} objects may be achieved with calls stated as
    follows: \code{new("saturateBF", bf.pv, nuLL, altLL, noInfo, nsites)}.
}

\section{Slots}{
    \describe{
        \item{\code{bf.pv}:}{numeric vector of the Bayes factors estimated
            independently for each molecular data site.}
        \item{\code{nuLL}:}{log-likelihood recorded for each data site based
            on the null parameters of the corresponding saturation analyses.}
        \item{\code{altLL}:}{estimates of the log-likelihood obtained for each
            data site, using the parameters of the alternative hypothesis.}
        \item{\code{noInfo}:}{site indices of the analysed protein data sites
            that are observed to be invariant.}
        \item{\code{nsites}:}{number of base sites that make up the original
            molecular data.}
    }
}

\section{Methods}{
    \describe{
        \item{show:}{\code{signature(object="saturateBF")} self-explanatory
            output highlights.}
        \item{summary:}{\code{signature(object="saturateBF")} interesting
            Bayesian-informed saturation statistics in properly labelled
            data frame.}
        \item{nonvaries:}{\code{signature(cNS="saturateBF")} vector that
            contains indices of the sites in the analysed protein data that
            were invariant.}
        \item{sitecount:}{\code{signature(gent="saturateBF")} number of base
            locations that makeup the molecular data used for the saturation
            analyses.}
        \item{BFs:}{\code{signature(satube="saturateBF")} Bayes factor
            estimates obtained for each data site.}
        \item{LogL0:}{\code{signature(satube="saturateBF")} estimates of the
            log-likelihood computed using the parameters of the null
            hypothesis, for each base site separately.}
        \item{LogL1:}{\code{signature(satube="saturateBF")} log-likelihood
            estimates for each data site as obtained from application of the
            parameters of the alternative hypothesis.}
    }
}

\section{Details}{
    The \code{\link{saturatePV}} class is a parent class, where the
    input of the \code{bf.pv} slots has a different interpretation.
}

\author{
    Hassan Sadiq
}

\seealso{
    The parent \code{\link{saturatePV}} class that is applicable to
    frequentist-based outputs of site-wise molecular saturation analyses.
    It may also be useful to see the \code{\link{seqSaturation}} function
    that produces objects of this class.
}

\examples{
    geneFile <- queeemsExtdata("II.fasta")
    examine <- seqSaturation(geneFile, 1)
    examine
}

\keyword{classes}

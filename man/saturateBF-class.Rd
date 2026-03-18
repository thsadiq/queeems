\name{saturateBF-class}
\docType{class}
\alias{saturateBF}
\alias{saturateBF-class}
\alias{BFs}
\alias{LogL0}
\alias{LogL1}
\alias{mainPi}
\alias{BFthreshold}
\alias{BFs,saturateBF-method}
\alias{show,saturateBF-method}
\alias{LogL0,saturateBF-method}
\alias{LogL1,saturateBF-method}
\alias{mainPi,saturateBF-method}
\alias{summary,saturateBF-method}
\alias{sitecount,saturateBF-method}
\alias{nonvaries,saturateBF-method}
\alias{BFthreshold,saturateBF-method}

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
        \item{\code{bfcutoff}:}{Bayes factor threshold for determining
            significance.}
        \item{\code{sspi}:}{tolerated proportion of saturated sites to use
            for overall hypotheses statements.}
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
        \item{mainPi:}{\code{signature(satube="saturateBF")} a positive
            fraction that specifies the value to use when conducting
            hypotheses tests for overall saturation.}
        \item{BFthreshold:}{\code{signature(satube="saturateBF")} Bayes
            factor threshold value that determines when to infer saturation.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    The \code{\link{seqSaturation}} function
    that produces objects of this class.
}

\examples{
    geneFile <- queeemsExtdata("II.fasta")
    examine <- seqSaturation(geneFile, "A", 5, 0.4, 5)
    examine
}

\keyword{classes}

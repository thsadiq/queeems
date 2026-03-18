\name{siteindices-class}
\docType{class}
\alias{siteindices}
\alias{siteindices-class}
\alias{sitentropies}
\alias{show,siteindices-method}
\alias{summary,siteindices-method}
\alias{renyiA,siteindices-method}
\alias{nonvaries,siteindices-method}
\alias{informula,siteindices-method}
\alias{sitentropies,siteindices-method}

\title{Molecular Information Entropy Processor}

\description{
    Allow for neat and efficient handling of information entropy indices
    outputs generated site-wise from molecular data.
}

\section{Objects from the Class}{
    To create a \code{siteindices} object require executing a call of the
    form \code{new("siteindices", eindices, noinfo, meth, alphaR)}.
}

\section{Slots}{
    \describe{
        \item{\code{eindices}:}{vector of site-wise entropy estimates.}
        \item{\code{noinfo}:}{numeric vector that contains positions of
            the molecular data that produced \code{eindices} which were
            non-informative.}
        \item{\code{meth}:}{character entry that informs about the method
            used to generate \code{eindices} from the original protein data.}
        \item{\code{alphaR}:}{parameter for Renyi entropy calculations. Only
            relevant if \code{meth = "renyi"}.}
    }
}

\section{Methods}{
    \describe{
        \item{show:}{\code{signature(object="siteindices")} output highlights.}
        \item{summary:}{\code{signature(object="siteindices")} interesting
            output statistics.}
        \item{nonvaries:}{\code{signature(cNS="siteindices")} positions of the
            processed molecular data that are not informative.}
        \item{renyiA:}{\code{signature(sent="siteindices")} value of the
            alpha parameter used for Renyi entropy calculations.}
        \item{informula:}{\code{signature(sent="siteindices")} should
            return the name of the method used to generate the entropy
            estimates from the original molecular data.}
        \item{sitentropies:}{\code{signature(sent="siteindices")} values of
            entropy estimates obtained for each site of a genetic sequence.}
    }
}

\author{
    Hassan Sadiq
}

\section{Details}{
    An approximate 95\% confidence interval estimated using the
    Bienayme-Tchebcshef inequality is returned as part of the outputs
    of the \code{siteindices} class. The expression for the said confidence
    interval calculation is as follows.
    \deqn{
        \bar{H} \pm k \times S_{H}^{} / \sqrt{n},
    }
    where \eqn{\bar{H}} and \eqn{S_{H}^{}} are the average and the standard
    deviation of the \eqn{n} entropy estimates supplied as the input for
    the \code{eindices} slot. The parameter \eqn{k=\sqrt{1/0.05}}.
}

\seealso{
    \itemize{
        \item Variant of the \code{nonvaries} method as applied within
            the \code{\link{baseSummary}} and the \code{\link{geneindex}}
            classes.
        \item The \code{\link{entropyindex}} function is a good reference
            for discussions about the different entropy calculation
            techniques accessible from the \code{\link{queeems}} package.
        \item Functions that return outputs of the described
            \code{siteindices-class} include \code{\link{codondifferindex}}
            and \code{\link{molentropy}}.
    }
}

\references{
    Feller, W. (1968). An Introduction to Probability Theory and its
    Applications: Vol. 1, \emph{Wiley New York}.
}

\examples{
    proteinFile <- queeemsExtdata("VertCOI.fasta")
    testEG <- codondifferindex(proteinFile, 1, "shannon")
    sitentropies(testEG)[seq(1,10)]
    nonvaries(testEG)[seq(1,10)]
    summary(testEG)
    testEG
}

\keyword{classes}

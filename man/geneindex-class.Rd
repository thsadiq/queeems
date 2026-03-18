\name{geneindex-class}
\docType{class}
\alias{geneindex}
\alias{geneindex-class}
\alias{renyiA}
\alias{useSyn}
\alias{gentropy}
\alias{informula}
\alias{sitecount}
\alias{show,geneindex-method}
\alias{renyiA,geneindex-method}
\alias{useSyn,geneindex-method}
\alias{gentropy,geneindex-method}
\alias{nonvaries,geneindex-method}
\alias{informula,geneindex-method}
\alias{sitecount,geneindex-method}

\title{Protein Information Entropy Handling}

\description{
    Facilitates succinct handling of the global information entropy
    index estimated for an alignment of homologous protein sequences.
}

\section{Objects from the Class}{
    A \code{geneindex} object may be created with a call of the form
    \code{new("geneindex", eindex, noinfo, meth, alphaR, nbases, synonym)}.
}

\section{Slots}{
    \describe{
        \item{\code{eindex}:}{single numeric entropy estimates obtained
            from genetic data.}
        \item{\code{noinfo}:}{numeric vector that contains non-informative
            positions of the molecular data that produced \code{eindex}.}
        \item{\code{meth}:}{character entry that informs about the method
            used to generate \code{eindex} from the original protein data.}
        \item{\code{alphaR}:}{parameter for Renyi entropy calculations. Only
            relevant if \code{meth = "renyi"} \code{[default = NA]}.}
        \item{\code{nbases}:}{number of base positions that make up the
            protein data used to estimate \code{eindex}.}
        \item{\code{synonym}:}{\code{default = NA} and only relevant if the
            counts used to estimate \code{eindex} are with respect to
            [non-]synonymous counts.}
    }
}

\section{Methods}{
    \describe{
        \item{show:}{\code{signature(object="geneindex")} output highlights.}
        \item{renyiA:}{\code{signature(sent="geneindex")} value of the
            alpha parameter used for Renyi entropy calculations.}
        \item{useSyn:}{\code{signature(gent="geneindex")} answers the
            question: were synonymous codon counts used?.}
        \item{gentropy:}{\code{signature(gent="geneindex")} value of the
            entropy estimate obtained for the genetic sequence.}
        \item{nonvaries:}{\code{signature(cNS="geneindex")} positions of the
            processed molecular data that are not informative.}
        \item{informula:}{\code{signature(sent="geneindex")} should
            return the name of the method used to generate the entropy
            estimate from the original molecular data.}
        \item{sitecount:}{\code{signature(gent="geneindex")} number of
            base positions that constitute the analysed protein data.}
    }
}

\author{
    Hassan Sadiq
}

\section{Details}{
    The entry of the \code{synonym} slot is expected to be logical, if
    specified. \code{TRUE} is set if the codon counts used to estimate
    the \code{eindex} entry is with respect to synonymous counts and
    should be \code{FALSE} if non-synonymous counts are used.
}

\seealso{
    \itemize{
        \item Variant of the \code{nonvaries} method is applied within
            the \code{\link{baseSummary}} and the \code{\link{siteindices}}
            classes.
        \item The \code{\link{entropyindex}} function is a good reference
            for discussions about the different entropy calculation
            techniques accessible from the \code{\link{queeems}} package.
        \item Functions that return outputs of the described
            \code{geneindex-class} include \code{\link{cncsentropy}}.
    }
}

\examples{
    geneFile <- queeemsExtdata("PI.fasta")
    testing <- cncsentropy(geneFile, FALSE, 1, "renyi", 0.04)
    print( nonvaries(testing))
    print( class(testing))
    print(testing)
}

\keyword{classes}

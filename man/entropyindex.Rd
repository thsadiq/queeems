\name{entropyindex}
\alias{entropyindex}
\alias{si.shannon}
\alias{si.renyi}

\title{Obtain Entropy Indices}

\description{
    Obtain indices that are useful for quantifying the degree of diversity
    contained in observed molecular prevalence count data.
}

\usage{
    entropyindex(countvector, nleaf, fml="shannon", ralpha=0.1)
}

\arguments{
    \item{countvector}{integer vector of base frequencies.}
    \item{nleaf}{number of the extant branches of the tree that generated
        the data reduced to \code{countvector}.}
    \item{fml}{character input that can be one of \code{"shannon"} or
        \code{"renyi"}.}
    \item{ralpha}{a fraction between 0 and 1 that is needed for Renyi
        index calculation.}
}

\details{
    Molecular entropy estimation requires knowledge of the size of the leaves
    on the phylogeny that generated the data. This is expected to be entered
    into the function as \code{nleaf}. Let \eqn{\bold{x}} denote the vector
    of base frequencies \code{countvector} and let \code{n = nleaf}. If
    \code{fml} is set as \code{"shannon"}, the diversity index is computed as:
    \deqn{
        H = -\sum_{i=1}^{m} p_{i} \times \log p_{i},
    }
    where \code{m} \eqn{\leq} \code{n} is the length of \eqn{\bold{x}} after
    all zeros have been excluded and \eqn{\code{p}_\code{i}} \code{=}
    \eqn{\code{x}_\code{i} / \code{n}}. When \code{fml} is set as
    \code{"renyi"}, then
    \deqn{
        H = \frac{1}{1-\alpha} \log \sum_{i=1}^{m} p_{i}^{\alpha},
    }
    where \eqn{\alpha} \code{= ralpha}.
}

\value{
    A numeric estimate of the entropy index obtained for the input data.
}

\author{
    Hassan Sadiq
}

\seealso{
    Scaled multinomial likelihood estimator \code{\link{mnomLogl}} which
    returns a scaled entropy. The \code{\link{seqSaturation}} relies heavily
    on this function, while \code{\link{baseFrequency}}, \code{\link{CnCs}}
    and \code{\link{codonDissimilarity}} are useful for transforming genetic
    sequences into counts needed for the \code{countvector} input. Interested
    users might consult the external \code{\link{vegan::diversity}} function
    for alternative entropy estimation methods.
}

\references{
    Renyi, A. (1961). On Measures of Information and Entropy,
    \emph{Proceedings of the Fourth Berkeley Symposium on Mathematics,
    Statistics and Probability 1960} 547-561.

    Shannon, C. E. (1948). A Mathematical Theory of Communication, \emph{The
    Bell System Technical Journal} \bold{27}(3): 379-423.
}

\examples{
    freqdata <- sample(seq(0,6), 20, replace=TRUE)
    tested <- entropyindex(freqdata, 20)
    tested
}

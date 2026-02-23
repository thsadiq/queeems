\name{mnomLogl}
\alias{mnomLogl}

\title{Approximate Log-Likelihood of a Multinomial Variable}

\description{
    Obtain an approximate estimate of the log-likelihood for a random
    variable that is believed to have a multinomial distribution.
}

\usage{
    mnomLogl(countvector, pivector)
}

\arguments{
    \item{countvector}{vector of non-negative observed frequencies.}
    \item{pivector}{theoretical population proportions that generated
        the input \code{countvector}.}
}

\details{
    Consider a discrere random variable
    \deqn{
        \mathbf{X} = \left(X_{1},\ldots,X_{k}\right) \sim
        \texttt{Multinomial}\left(n,\,\pi_{1},\ldots,\pi_{k}\right),
    }
    where, \eqn{n = X_{1} + \ldots + X_{k}} and \eqn{\sum_{i=1}^{k}
    \pi_{i} = 1}. Then,
    \deqn{
        P\left(\mathbf{x}=\left(x_{1},\ldots,x_{k}\right)\,|\,
            \boldsymbol{\pi} \right) = \frac{n!}{x_{1}!,\ldots,x_{k}!}\,
            \pi_{1}^{x_{1}} \ldots,\pi_{k}^{x_{k}}
    }
    such that the log-likelihood expression may be simplified as follows.
    \deqn{
        l\left(\boldsymbol{\pi}\,|\,\mathbf{x}\right) = \log \left( n! \right)
        - \sum_{i=1}^{k} \log \left( x_{i}! \right) + \sum_{i=1}^{k} x_{i}
        \log \left( \pi_{i} \right) \; \propto \; \sum_{i=1}^{k} x_{i}
        \log \left( \pi_{i} \right).
    }
    The expression can be shown to produce a maximum likelihod estimate of
    \eqn{\hat{\pi}_{i} = x_{i} / n} by using a Lagrange multiplier that
    imposes the constraint that \eqn{\sum_{i=1}^{k} \pi_{i} = 1}.. 
}

\value{
   The scaled likelihood value estimated from the final expression presented
   in \bold{Details}. The output should be a negative value.
}

\author{
    Hassan Sadiq
}

\seealso{
    Some other information retrieval functions within the \code{\link{queeems}}
    package that are straightforwardly applicable to molecular protein data.
    These include, \code{\link{dnaEntropy}}, \code{\link{cncsEntropy}} and
    \code{\link{seqSaturation}}. Functions such as \code{\link{baseFrequency}},
    \code{\link{CnCs}} and \code{\link{codonDissimilarity}}, that are needed
    to obtain the counts required to successfully execute the described
    \code{mnomLogl}, from a genetic sequences may also be consulted.
}

\examples{
    arbitraryCounts <- sample(40, 20)
    tempPis <- runif(20)
    arbitraryPis <- tempPis / sum(tempPis)
    tests <- mnomLogl(arbitraryCounts, arbitraryPis)
    print(tests)
}

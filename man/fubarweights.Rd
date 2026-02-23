\name{fubarweights}
\alias{fubarweights}

\title{FUBAR-like Weights}

\description{
    Generate unevenly discretised weights
}

\usage{
    fubarweights(nweights, wless1, weightmax, addzero)
}

\arguments{
    \item{nweights}{number of weights.}
    \item{wless1}{proportion of weights between the (0, 1) interval.}
    \item{weightmax}{maximum weight value.}
    \item{addzero}{logical input. \code{TRUE} if minimum weight value should
        be zero and \code{FALSE} otherwise.}
}

\details{
    For \eqn{n > 3} unevenly discretised weights that are all less than
    \eqn{w_{max}} such that \eqn{\phi} of them lie within the (0, 1)
    interval, the vector of weights is obtained as the following set of
    values.
    \deqn{
        W = \{w_{L}^{}, w_{U}^{}\},
    }
    such that,
    \deqn{
        w_{L}^{} = \{\frac{I^{(0)}}{\kappa}, \frac{I^{(0)}+1}{\kappa},
        \ldots, \frac{\kappa}{\kappa}\},
    }
    \deqn{
        w_{U}^{} = \{1+(\delta)^3, 1+(2\delta)^3, \ldots,
        1+((n-\kappa+I^{(0)}-1)\delta)^3\},
    }
    where
    \deqn{
        \kappa = \lfloor \phi n \rfloor,
    }
    \deqn{
        \delta = \frac{\sqrt[3]{w_{max}-1}}{n-\kappa+I^{(0)}-1},
    }
    Further, \eqn{I^{(0)}} equals zero if \code{addzero = TRUE} and it equals
    one otherwise, \eqn{n} is the input \code{weightmax} and \eqn{\phi} is
    the specified \code{wless1}.
}

\value{
    A numeric vector of weights.
}

\author{
    Hassan Sadiq
}

\references{
    Murrell, B., Moola, S., Mabona, A., Weighill, T., Sheward, D.,
    Kosakovsky Pond, S. L. and Scheffler, K. (2013) FUBAR: A Fast,
    Unconstrained Bayesian AppRoximation for Inferring Selection,
    \emph{Mol. Biol. Evol.} \bold{30}(5): 1196-1205.
}

\examples{
    wVector1 <- fubarweights(20, 0.70, 50, FALSE)
    wVector2 <- fubarweights(20, 0.70, 50, TRUE)
    print(wVector1)
    print(wVector2)
}

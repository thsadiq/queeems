\name{softmax}
\alias{softmax}

\title{Softmax Transformation}

\description{
    Transform weights to frequencies using the softmax identity
}

\usage{
    softmax(weights, maxadj)
}

\arguments{
    \item{weights}{vector of weights.}
    \item{maxadj}{logical indicator of whether weights need adjustment.}
}

\details{
    The softmax identity is popular in machine learning applications. It is
    useful for transforming a vector of weights into frequencies. A modified
    version that minimises the chances of overflow is applied here. For a
    given vector (\eqn{x} = \code{weights}) that contains \eqn{k} real number
    weights, the \code{i}th frequency (where, \code{i = 1, 2, ..., k}) is
    obtained as follows.
    \deqn{
        \pi_{i}^{} =
            \frac{e^{x_{i}^{}-c}_{}}{\sum_{j=1}^{k} e^{x_{j}^{}-c}_{}},
    }
    where \eqn{c} is set equal to zero if \code{maxadj = FALSE}. If
    \code{maxadj = TRUE} then, \eqn{c=\max_{j} x_{j}^{}}.
}

\value{
    A numeric vector of frequencies that sum to one.
}

\author{
    Hassan Sadiq
}

\references{
    Efron, B. and Hastie, T. (2016). \emph{Computer Age Statistical
    Inference. Algorithms, Evidence and Data Science}. Cambridge, UK:
    Cambridge University Press.
}

\examples{
    wData <- sample(seq(1,6), 20, replace=TRUE)
    test1 <- softmax(wData, FALSE)
    test2 <- softmax(wData, TRUE)
    print(test1)
    print(test2)
}

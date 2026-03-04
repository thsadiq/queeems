\name{aboutQueeems}
\alias{aboutQueeems}

\title{References Relevant to the \code{queeems} Package}

\description{
    Obtain citations of research outputs that may assist with understanding
    and extending the applications of the \code{queeems} package.
}

\usage{
    aboutQueeems(cite=NULL)
}

\arguments{
    \item{cite}{A strictly positive integer.}
}

\details{
    Citations included are with respect to outputs that significantly use
    the \code{queeems} package and that the package author contributes to.
    The list will be extended as more outputs are published. When
    \code{cite = NULL} is set (default), all the available citations are
    listed. If the integer input for \code{cite} exceeds the number of
    available citations, an error message is returned.
}

\value{
    A text that contains citations associated with the \code{queeems}
    package, as explained in \bold{Details}.
}

\author{
    Hassan Sadiq
}

\examples{
    aboutQueeems()
}

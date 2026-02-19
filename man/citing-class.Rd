\name{citing-class}
\docType{class}
\alias{citing}
\alias{citing-class}
\alias{show, citing-method}

\title{Useful Package Citations}

\description{
    Get details of citations that may further assist with understanding of
    package utility.
}

\section{Objects from the Class}{
    Objects of this (\code{citing}) class can be created by calls of the
    form \code{new("citing", citeText=...)}.
}

\section{Slots}{
    \describe{
        \item{\code{citeText}:}{printable text that may be pasted directly
            as references into scientific write-ups.}
    }
}

\section{Methods}{
    \describe{
        \item{show:}{\code{signature(object = "citing")} prints the input
            text.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    The \code{\link{aboutQueeems}} function that relies on the described class.
}

\examples{
    test <- aboutQueeems()
    is(test, "citing")
    test
}


\keyword{classes}

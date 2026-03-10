\name{queeems}
\alias{queeems}
\alias{queeems-package}
\docType{package}

\title{Quantify the Extent of Evolutionary Evidence in Molecular Sequences}

\description{
    Saturation indices are often used, in molecular evolutionary studies, to
    determine whether the degree of genetic variation captured in a set of
    genetic sequences is likely to produce reliable inferences. Existing
    measures are largely designed to cater for analyses about phylogeny
    reconstruction. As a result, a single metric is usually obtained for
    an alignment of homologous sequences. When saturation is detected, it
    is common practice to exclude sequence positions believed to
    significantly influence the erosion of genetic variation. The suite of
    functions presented here are designed to facilitate a more informed
    exclusion technique by providing Bayesian saturation metric for each
    sequence position.
}

\details{
    Bayes factors that measure the extent to which the recorded genetic data
    support an alternative hypothesis of saturation is estimated for each
    sequence position, in terms of codon bases. This approach avails an
    exciting opportunity for saturation indices to find use in other
    molecular evolutionary topics, including analyses of natural selection.
    Other independently useful functions are also part of the package, these
    include various molecular entropy quantifying functions.
}

\author{
    Hassan Sadiq
}

\references{
    Kass, R. E. and Raftery, A. E. (1995). Bayes Factors, \emph{Journal of
    American Statistical Association} \bold{90}(430): 773-795.

    Renyi, A. (1961). On Measures of Information and Entropy,
    \emph{Proceedings of the Fourth Berkeley Symposium on Mathematics,
    Statistics and Probability 1960} 547-561.

    Shannon, C. E. (1948). A Mathematical Theory of Communication, \emph{The
    Bell System Technical Journal} \bold{27}(3): 379-423.

    Xia, X., Xie, Z., Salemi, M., Chen, L. and Wang, Y. (2003). An Index
    of Substitution Saturation and its Application, \emph{Molecular
    Phylogenetics and Evolution} \bold{26}(1): 1-7.
}

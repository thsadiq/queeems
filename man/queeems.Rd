\name{queeems}
\alias{queeems}
\alias{queeems-package}
\docType{package}

\title{Quantify the Extent of Evolutionary Evidence in Molecular Sequences}

\description{
    Saturation indices are often used, in molecular evolutionary studies, to
    determine whether the degree of genetic variation captured in a set of
    protein sequences is likely to produce reliable inferences. Existing
    measures are largely designed to cater for analyses about phylogeny
    reconstruction. As a result, a single metric is usually obtained for
    an alignment of homologous sequences. When saturation is detected, it
    is common practice to exclude sequence positions believed to
    significantly influence the erosion of genetic variation. The suite of
    functions presented here are designed to facilitate a more informed
    exclusion technique by providing saturation metric for each sequence
    position.
}

\details{
    Bayes factors that measure the extent to which the recorded protein data
    support an alternative hypothesis of no saturation is estimated for
    each sequence position, in terms of DNA, codon or amino acid bases. This
    approach avails an exciting opportunity for saturation indices to find
    use in other molecular evolutionary topics, including analyses of
    natural selection.
}

\author{
    Hassan Sadiq
}

\references{
    Kass, R. E. and Raftery, A. E. (1995). Bayes Factors, \emph{Journal of
    American Statistical Association} \bold{90}(430): 773-795.

    Xia, X., Xie, Z., Salemi, M., Chen, L. and Wang, Y. (2003). An Index
    of Substitution Saturation and its Application, \emph{Molecular
    Phylogenetics and Evolution} \bold{26}(1): 1-7.

    Sadiq, H. (in progress). queeems: Quantify the Extent of Evolutionary
    Evidence in Molecular Sequences, \emph{R Package}.

    Sadiq, H. (in progress). Saturation Threshold of Codon Evolutionary
    Models, \emph{Preprint}.
}

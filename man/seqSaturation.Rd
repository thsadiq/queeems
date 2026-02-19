\name{seqSaturation}
\alias{seqSaturation}

\title{Estimate Saturation Index For Protein Data}

\description{
    Estimate saturation index for molecular sequences
    with respect to DNA, codon or amino acid bases.
}

\usage{
    seqSaturation(fastafile, approach, gammaPset=c(l=0.1,u=10,mw=3))
}

\arguments{
    \item{fastafile}{a character input of the full path to a \code{.fasta}
        file that contains the protein sequence(s) to be analysed.}
    \item{basename}{one of `\code{dna}`, `\code{codon}` or `\code{aa}` that
        indicates how to decipher the bases of the input sequence.}
    \item{approach}{an integer that indicates the preferred saturation
        assessment method. It can be any one of \code{{1,2,3}}.}
}

\details{
    Three molecular saturation verification techniques are accessible
    through the \code{seqSaturation} function. The logic behind the
    techniques are the following.
    \itemize{
        \item \bold{Technique 1 (\code{approach = 1}):} codons are
            combinations of nucleotide bases. This attribute enable codons
            to be particularly useful for more extensive evolutionary
            analyses, including the popular investigation of natural
            selection as a ratio of non-synonymous to synonymous codon
            substitution rates. In addition, due to the smaller set of
            nucleotide compared to codon bases, the former are expected
            to lose evolutionary information faster. The substitution
            saturation test built upon this reasoning culminates in the
            following hypotheses statements.
        \itemize{
            \item \emph{Null hypothesis (\code{H0})}: Retrievable
                evolutionary information is similar for both nucleotide and
                codon data. That is, sequence is saturated.
            \item \emph{Alternative hypothesis (\code{H1})}: Evolutionary
                information contained in codon data is \bold{more} than the
                information retrievable from associated nucleotide data.
                That is, sequence is \bold{not} saturated.
        }
        These hypotheses are verified using Bayes factors generated from
        the following expression.
        \deqn{
            A = B.
        }

        \item \bold{Technique 2 (\code{approach = 2}):} a significant number
            of existing codon models of evolution tolerate only one
            nucleotide mismatch between codon pairs. The assumption is
            that not enough time has passed for more than a single
            substitution event to have occurred. The following hypotheses
            statements are consequence of this widely adopted assumption.
        \itemize{
            \item \emph{Null hypothesis (\code{H0})}: Observed base
                frequencies distribution is better supported by a model that
                permits both two and three nucleotide mismatches per codon
                pair. That is, sequence is saturated.
            \item \emph{Alternative hypothesis (\code{H1})}: The frequencies
                of the codon bases that make up the protein data are better
                described by a model that supports only a single nucleotide
                mismatch between codon pairs. That is, sequence is \bold{not}
                saturated.
        }
        The Bayes factor for verifying these statements can be expressed
        as follows.
        \deqn{
            C = D.
        }

        \item \bold{Technique 3 (\code{approach = 3}):} this is an extension
            of \code{Technique 2} that allows for the hypotheses to be nested
            such that a likelihood ratio test becomes applicable with the
            following updated hypothesis statements.
        \itemize{
            \item \emph{Null hypothesis (\code{H0})}: Including information
                about one, two and three nucleotide mismatches per codon
                pair better capture the evolutionary evidence that have
                accrued in the protein data. That is, sequence is saturated.
            \item \emph{Alternative hypothesis (\code{H1})}: Considering only
                one nucleotide mismatch per codon pair captures as much of
                the evidence that has accrued in the protein data. That is,
                sequence is \bold{not} saturated.
        }
        These hypotheses statements are assessed with the following
        likelihood ratio expression.
        \deqn{
            E = F.
        }
    }
    All these techniques analysed more thoroughly in an upcoming article
    that heavily uses functions from the \code{\link{queeems}} package.
}

\value{
   The output object returned is determined by the analyses technique.
   Bayesian outputs from \code{Technique 1} and \code{Technique 2} are
   returned as \code{\link{saturateBF}} objects. Outputs from the
   frequentist-based \code{Technique 3} are returned as
   \code{\link{saturatePV}} objects. Both objects contain details of
   the site-wise entropy estimates and other useful characteristic
   information. The respective function help files may be consulted
   for further details.
}

\author{
    Hassan Sadiq
}

\seealso{
    \itemize{
        \item Molecular entropy estimation functions:
            \itemize{
                \item \code{\link{entropyindex}}
                \item \code{\link{codondifferindex}}
                \item \code{\link{entropyindex}}
                \item \code{\link{molentropy}}
            }
        \item Multinomial likelihood generating function: 
            \itemize{
                \item \code{\link{mnomLogl}}
            }
        \item Functions useful for transforming molecular data to useful
            frequencies:
            \itemize{
                \item \code{\link{baseFrequency}}
                \item \code{\link{codonDissimilarity}}
            }
    }
}

\references{
    Xia, X., Xie, Z., Salemi, M., Chen, L. and Wang, Y. (2003). An Index
    of Substitution Saturation and its Application, \emph{Molecular
    Phylogenetics and Evolution} \bold{26}(1): 1-7.

    Renyi, A. (1961). On Measures of Information and Entropy,
    \emph{Proceedings of the Fourth Berkeley Symposium on Mathematics,
    Statistics and Probability 1960} 547-561.

    Shannon, C. E. (1948). A Mathematical Theory of Communication, \emph{The
    Bell System Technical Journal} \bold{27}(3): 379-423.
}

\examples{
    proteinFile <- queeemsExtdata("VertCOI.fasta")
    testEG <- seqSaturation(proteinFile, 1)
    summary(testEG)
    testEG
}

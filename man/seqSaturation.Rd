\name{seqSaturation}
\alias{seqSaturation}

\title{Estimate Saturation Index For Protein Data}

\description{
    Assess the degree of evolutionary saturation in molecular sequences.
}

\usage{
    seqSaturation(fastafile, approach, wsize, wless1,
        weightmax, bfcutoff=10, sspi=0.10)
}

\arguments{
    \item{fastafile}{a character input of the full path to a \code{.fasta}
        file that contains the protein data to be analysed.}
    \item{approach}{one of `\code{A}`, `\code{B}` or `\code{C}` that
        indicates the preferred saturation assessment method.}
    \item{wsize}{number of \eqn{\pi} parameter weight classes.}
    \item{wless1}{fraction of \eqn{\pi} weight classes to be less than one.}
    \item{weightmax}{maximum \eqn{\pi} weight value.}
    \item{bfcutoff}{Bayes factor cutoff.}
    \item{sspi}{tolerated proportion of saturated sites, overall.}
}

\details{
    At least, three molecular saturation verification techniques are planned
    to be incorporated into the package. However, only one is currently
    accessible through the \code{seqSaturation} function. The logic behind
    the technique is presented below.
    \itemize{
        \item \bold{Technique 1 (\code{approach = "A"}):}
        \itemize{
            \item Codons are combinations of nucleotide bases. Codons are
            therefore more appealing for rigorous evolutionary analyses,
            including the popular investigation of natural selection as a
            ratio of non-synonymous to synonymous codon substitution rates.
            In addition, due to the smaller set of nucleotide bases (= 4)
            compared to codon bases (= 64), the former are expected to lose
            evolutionary information faster. The substitution saturation test
            built upon this reasoning reduces to the following hypotheses
            statements.
            \item \emph{Null hypothesis (\code{H0})}: Evolutionary information
                contained in codon data is significantly \bold{more} than the
                information retrievable from associated nucleotide data. That
                is, sequence is \bold{not} saturated. (\bold{\code{Bayes
                Factor}} \eqn{\leq} \bold{\eqn{\code{BF}_\code{cutoff}}}).
            \item \emph{Alternative hypothesis (\code{H1})}: Evolutionary
                information retrievable in the codon data is significantly
                eroded or insufficient for evolutionary analyses, such that
                nucleotide models reveal more about the data. That is,
                sequence is saturated. (\bold{\code{Bayes Factor}} \eqn{>}
                \bold{\eqn{\code{BF}_\code{cutoff}}}).

            \eqn{\code{BF}_\code{cutoff}} is the input \code{bfcutoff}. To
            verify these hypotheses, it is assumed that the number of codon
            and nucleotide bases observed at each of the corresponding
            positions of protein data follows a multinomial distribution.
            Bayes factors are then estimated independently for each codon
            site, \code{i = (1, 2, ..., n)}, as the ratio of maximum
            posterior log-likelihoods as follows.
            \deqn{
                \code{BF}_{i} = \frac{LogL\big(X_{\star i}^{(k)} |
                        \boldsymbol{\pi}_{\star i}^{(k)}\big)}
                        {LogL\big(Y_{\star i} |
                        \boldsymbol{\pi}_{\star i}^{\code{cdn}}\big)}
                    \propto
                    \frac{\displaystyle{\sum_{k=1}^{3}} \bigg(
                    \max_{\boldsymbol{\pi}^{(k)}}
                    \displaystyle{\sum_{j=1}^{4}} x_{ji}^{(k)} \times
                    \log\left(\pi_{ji}^{(k)}\right) \bigg)}
                    {\displaystyle{\max_{\boldsymbol{\pi}^{\code{cdn}}}}
                    \displaystyle{\sum_{m=1}^{61}} y_{mi} \times
                    \log\left(\pi_{mi}^{\code{cdn}}\right) },
            }
            where,
            \itemize{
                \item \eqn{Y_{\star i}} = data at codon position \code{i},
                \item \eqn{\boldsymbol{\pi}_{\star i}^{\code{cdn}}} =
                    theoretical proportions of sense codons at site \code{i},
                \item \eqn{X_{\star i}^{(k)}} = nucleotide data in position
                    \code{k = (1, 2, 3)} at codon site \code{i},
                \item \eqn{\boldsymbol{\pi}_{\star i}^{(k)}} = theoretical
                    proportions of nucleotide bases at position \code{k} of
                    codon site \code{i},
                \item \eqn{y_{mi}} = number of sense codon \code{m} observed
                    at codon \code{i},
                \item \eqn{x_{ji}^{(k)}} = number of nucleotide base \code{j}
                    recorded at position \code{k} of codon site \code{i}.
            }

            Values of \eqn{\boldsymbol{\pi}} that the likelihoods are maximised
            over for each site are pre-determined based on a technique adopted
            from FUBAR. Weights of size \code{w = nweights} \eqn{\geq 4} within
            the \code{(0, weightmax)} interval are first generated such that
            \code{wless1} fraction of the weights are less than one. For more
            details about the weight generation process, see
            \code{\link{fubarweights}}. All possible, \eqn{z = (w!)/(w-4)!},
            four-weight combinations without repetitions are then sampled and
            converted to \code{z} vectors of nucleotide frequencies
            \eqn{\left( \boldsymbol{\pi}^{(k)} \right)} using the
            \code{\link{softmax}} transformation. The \eqn{z} codon frequencies
            \eqn{\left( \boldsymbol{\pi}^{\code{cdn}} \right)} considered at
            each site are obtained as the product of the frequencies of the
            associated nucleotide triplets as it is conventionally implemented
            in the \code{F3x4} substitution models.

            The log-likelihood expressions are similar to the popular Shannon
            information entropy expression. The similarity justifies why it
            is not inappropriate to have constructed the hypotheses statements
            in terms of the information that can be retrieved. Note that the
            log-lilihood estimates will always be negative for both the
            numerator and denominator of the Bayes factor formula, because
            the \code{x, y} counts are strictly non-negative and the
            theoretical proportions (\eqn{\pi}) are always fractions.
            Consequently, to avoid spurious outcomes, the log-likelihood
            estimates were transformed into approximate posterior
            probabilities using the \code{\link{softmax}} identity before
            maximising.

            A single Bayes factor is also estimated to assess the overall
            extent of saturation of the protein sequence data (Y). Separate
            hypotheses statements about the true proportion of saturated
            sites (\eqn{\rho}), were necessary for this purpose. These
            hypotheses may be expressed as follows.
            \deqn{
                \code{H}_{0}:  \rho \leq \rho_{0}^{} \\
                \code{H}_{1}:  \rho > \rho_{0}^{},
            }

            where \eqn{\rho} is the parameter of the \code{Binomial}(\code{n},
            \eqn{\rho}) distribution, \code{n} = the number of codon sites and
            \eqn{\rho_{0}} is the input \code{sspi} (default \code{0.10}).
            Similar to the approach taken for the site-wise inferences,
            potential \eqn{\rho \in (0,1)} values are obtained using the
            FUBAR discretisation method. A Bayes factor is then estimated
            as the ratio of the maximum posterior probabilities as follows.
            \deqn{
                \code{BF}_{\code{overall}} = \frac{L \left(Y \big| \rho >
                    \rho_{0}^{} \right)}{L \left(Y \big| \rho \leq
                    \rho_{0}^{}\right)}.
            }

        }
    }
}

\value{
    The output object returned is determined by the analyses technique.
    Additional techniques should be incorporated soon. Bayesian outputs
    from \code{Technique 1} is returned as a \code{\link{saturateBF}} object.
    The object contains details of the site-wise likelihood estimates and
    other useful characteristic information. The respective object help
    file may be consulted for further details.
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
        \item Functions useful for transforming molecular data into the
            \code{x, y} counts used in the equations above:
            \itemize{
                \item \code{\link{baseFrequency}}
                \item \code{\link{codonDissimilarity}}
            }
        \item Weight generating function.
            \itemize{
                \item \code{\link{fubarweights}}
            }
        \item Function that converts nucleotide to codon frequencies.
            \itemize{
                \item \code{\link{n2cFreqs}}
            }
        \item Softmax identity function.
            \itemize{
                \item \code{\link{softmax}}
            }
    }
}

\references{
    Benjamini, Y. and Hochberg, Y. (1995). Controlling the False Discovery
    Rate: A Practical and Powerful Approach to Multiple Testing,
    \emph{Journal of the Royal Statistical Society: Series B}. \bold{57}(1):
    289-300.

    Efron, B. and Hastie, T. (2016). \emph{Computer Age Statistical
    Inference. Algorithms, Evidence and Data Science}. Cambridge, UK:
    Cambridge University Press.

    Goldman, N. and Yang, Z. (1994) A Codon-Based Model of Nucleotide
    Substitution for Protein-Coding DNA Sequences, \emph{Molecular
    Biology and Evolution} \bold{11}(5): 725-736.

    Kass, R. E. and Raftery, A. E. (1995) Bayes Factors, \emph{Journal
    of the American Statistical Association} \bold{90}(430): 773-795.

    Murrell, B., Moola, S., Mabona, A., Weighill, T., Sheward, D.,
    Kosakovsky Pond, S. L. and Scheffler, K. (2013) FUBAR: A Fast,
    Unconstrained Bayesian AppRoximation for Inferring Selection,
    \emph{Molecular Biology and Evolution} \bold{30}(5): 1196-1205.

    Renyi, A. (1961). On Measures of Information and Entropy,
    \emph{Proceedings of the Fourth Berkeley Symposium on Mathematics,
    Statistics and Probability 1960} 547-561.

    Shannon, C. E. (1948). A Mathematical Theory of Communication, \emph{The
    Bell System Technical Journal} \bold{27}(3): 379-423.

    Xia, X., Xie, Z., Salemi, M., Chen, L. and Wang, Y. (2003). An Index
    of Substitution Saturation and its Application, \emph{Molecular
    Phylogenetics and Evolution} \bold{26}(1): 1-7.
}

\examples{
    proteinFile <- queeemsExtdata("VertCOI.fasta")
    testEG <- seqSaturation(proteinFile, "A", 5, 0.4, 5)
    summary(testEG)
    testEG
}

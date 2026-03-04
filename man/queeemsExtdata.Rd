\name{queeemsExtdata}
\alias{queeemsExtdata}
\alias{II.fasta}
\alias{PI.fasta}
\alias{RTI.fasta}
\alias{VertCOI.fasta}

\title{Paths to Example Protein Sequence Files}

\description{
    Obtain the paths to example protein sequences that are available,
    as part of, and are used in some of the examples in, the \code{queeems}
    package. 
}

\usage{
    queeemsExtdata(filename=NULL)
}

\arguments{
    \item{filename}{A string of the name of one of the \code{.fasta} files
        packaged with \code{queeems}.}
}

\details{
    The external data sets accessible with this function include the,
    \code{II.fasta}, \code{PI.fasta} and \code{RTI.fasta} DNA files as
    well as their amino acid translates. That is, \code{II_AA.fasta},
    \code{PI_AA.fasta} and \code{RTI_AA.fasta}. The FASTA format
    files all contain empirical HIV-1 protein sequences that were retrieved
    from Murrell et al. (2012). Other files accessible from the folder
    include the \code{VertCOI.fasta}: mitochondrial COI sequences from
    eight vertebrates (Xia, 2018). The \code{filename} input may be
    specified as any of the listed filnames or left as the default
    \code{NULL} value. 
    
    The reverse transcriptase isolates [\code{RTI.fasta}] contain 205 codon
    sites and 476 sequences. Following Murrell et al. (2012), this corresponds
    to a subset [codon sites 41 to 245] of the original (Rhee et al, 2003,
    de Oliveira et al,, 2010) data. The protease isolate [\code{PI.fasta}]
    file contains 99 codon sites and 122 protein sequences. The integrase
    isolates [\code{II.fasta}] comprise 295 sequences and 288 codon sites.
    All the cited HIV-1 data sets contain pre- and post-treatment sequences.
}

\value{
    A string that corresponds to the full path of the data saved with the
    specified \code{filename} input. If \code{filename} is left as its default
    \code{NULL} value, a list of all the names of all the sequence files that
    are included with the \code{queeems} package are listed.
}

\author{
    Hassan Sadiq
}

\seealso{
    See applications of this function as part of \code{\link{bstringCodons}}.
}

\references{
    Rhee, S.-Y., Gonzales, M. J., Kantor, R., et al. (2003). Human
    Immunodeficiency Virus Reverse Transcriptase and Protease Sequence
    Database, \emph{Nucleic Acids Research}
    \bold{31}(1):298-303.

    de Oliveira, T., Shafer, R. W., Seebregts, C., et al. (2010). Public
    Database for HIV Drug Resistance in Southern Africa, \emph{Nature}
    \bold{464}(7289):673.

    Murrell, B., de Oliveira, T., Seebregts, C., et al. (2012). Modeling HIV-1
    Drug Resistance as Episodic Directional Selection, \emph{PLoS Computational
    Biology} \bold{8}(5):e1002507.

    Xia, X. (2018). DAMBE7: New and Improved Tools for Data Analysis in
    Molecular Biology and Evolution, \emph{Molecular Biology and Evolution}
    \bold{35}(6):1550-1552.
}

\examples{
    test0 <- queeemsExtdata()
    print(test0)

    test1 <- queeemsExtdata("II.fasta")
    print(test1)
}

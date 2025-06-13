---
title: 'scoup: Simulate Codon Sequences with Darwinian
        Selection Incorporated as an Ornstein-Uhlenbeck Process'
tags:
  - R
  - Molecular Evolution
  - Ornstein-Uhlenbeck
  - Phylogenetics
  - Population Genetics
  - Simulation
  - Statistics
authors:
  - name: Hassan Sadiq
    orcid: 0000-0003-0192-7134
    affiliation: "1, 2"
  - name: Darren P. Martin
    orcid: 0000-0002-8785-0870
    affiliation: 2
affiliations:
 - name: Department of Statistics and Actuarial Science,
         Stellenbosch University, South Africa
   index: 1
 - name: Institute of Infectious Diseases and Molecular Medicine,
         Division of Computational Biology,
         Department of Integrative Biomedical Sciences,
         University of Cape Town, South Africa
   index: 2
date: 13 June 2025
bibliography: paper.bib
---

# Summary

Genetic analyses of natural selection within and between populations have
increasingly developed along separate paths. The two important genres of
evolutionary biology (i.e. phylogenetics and population genetics) borne
from the split can only benefit from research that seeks to bridge
the gap. Simulation algorithms that combine fundamental concepts from
both genres are important to achieve such unifying objective. 
We introduce `scoup`, a codon sequence simulator that is implemented
in R and hosted on the Bioconductor platform. There is hardly any other
simulator dedicated to genetic sequence generation for natural selection
inference on the platform. Concepts from the Halpern-Bruno mutation-selection
model and the Ornstein-Uhlenbeck (OU) evolutionary algorithm were creatively
fused such that the end-product is a novelty with respect to computational
genetic simulation. Users are able to seamlessly adjust the model parameters
to mimic complex evolutionary procedures that may have been otherwise
infeasible. For example, it is possible to explicitly interrogate the
concepts of static and changing fitness landscapes with regards to Darwinian
natural selection in the context of codon sequences from multiple
populations.

# Statement of need

Statistical inference of the extent to which Darwinian natural selection has
impacted genetic data commands a healthy portion of the phylogenetic
literature [@jacques2023]. Validation of these largely codon-based models
relies heavily on simulated data. Given the ever increasing diversity of
natural selection inference models that exist [@yang2007; @hyphy2020],
there is a need for more sophisticated simulators to match the ever
expanding model complexities.

Bioconductor [@bioc2004] is a leading platform where peer-reviewed
bioinformatic software useful for biological data analyses are hosted. A
search of the entries on the platform, in Version 3.19 on 29 October 2024,
with keywords including, `codon`, `mutation`, `selection`, `simulate`, and
`simulation` returned a total of 72 unique packages out of the 2300 available.
None of the retrieved entries was dedicated to codon data simulation for
natural selection analyses. Thus, `scoup` is designed on the basis of the
mutation-selection (MutSel) framework [@halpern1998] as an overdue
contribution to the void. Software and/or packages for simulating genetic
sequences are also rare in the scientific literature [@gearty2024]. 

# Algorithm

`scoup` is further unique for at least three reasons. First, it incorporates
Darwinian natural selection into the MutSel model in terms of variability of
selection coefficients [@spielman2015]. Second, it directly utilises the
concept of fitness landscapes. Third, fitness landscape updates can be
executed in either a deterministic or a stochastic format. The stochastic
updates are implemented in terms of the more biologically amenable,
Ornstein-Uhlenbeck (OU) process [@bartoszek2017; @uhlenbeck1930]. A crude
summary of how substitution events are initiated and progresses in
`scoup` is presented in Figure \autoref{sfrrame}.

![\label{sfrrame}**Summarised `scoup` algorithm.**
  After each substitution event and until the input tree length
  ($\tau \in \mathbf{SEQ}$) is exhausted, the process returns to
  *STEP A*. $\sigma^{2}_{`n`}=$ variance of amino acid selection
  coefficients. $\sigma^{2}_{`s`}=$ variance of synonymous codon
  selection coefficients. $\Sigma^{2}_{}=$ OU asymptotic variance.
  $\theta=$ OU mean reversion rate. $\mathbf{SEQ}=$ sequence
  alignment details. $x_{\star}^{}=$ codon. $\mathbf{s}_{\star}^{}=$
  codon selection coefficient vector.](FIG0.pdf){ width=58% }

# Implementation

We simulated (see sample code in Figure \ref{pseudocode}) $20$ independent
sequence alignments made up of $1000$ codon sites and $8$ extant taxa for
each of the parameter combinations presented in Figure \ref{testfig}. The
phylogeny used was balanced and the length of its branches were $0.10$ each.
The stochastic OU framework was implemented and other function inputs were
left at their default values (Figure \ref{pseudocode}). Estimates of
$\mathrm{d}N/\mathrm{d}S$ were obtained following @spielman2015 and were
averaged over all selection coefficient updates at each site and across
the alignment. Inferences of $\omega$ were obtained with `CODEML` in
`PAML` [@yang2007].

![\label{pseudocode}**An example R code for simulating a codon sequence
  alignment with `scoup`**. Default values were left unchanged. `Line01`:
  OU adaptation parameters where, $\mu=0$, $\Sigma^{2}_{}=0.01$ and
  $\theta=0.01$. `Line02`: evolution model input where,
  $\mathbf{s} \sim \text{Gamma}(1,\sigma_{n}^{-1})$,
  $\sigma^{2}_{`n`}=10^{-5}_{}$, $\sigma^{2}_{`s`}=10^{-5}_{}$ and
  effective population size, $N_{\texttt{e}}^{}=1000$. `Line03`:
  sequence information where, site count is $250$, extant taxa count
  is $64$ and branch length is $0.1$.](FIG1.png){ width=50% }

![\label{testfig}**Demonstration of the accuracy of outputs from `scoup`
  in terms of the likelihood $\omega$ and the analytical dN/dS measures
  of natural selection.** The estimates of the selection measures were
  obtained homogeneously from each alignment generated for every combination
  of the stochastic landscape ($\Sigma^{2}_{}$ and $\theta$) and the
  Darwinian selection ($\sigma^{2}_{`n`}$ and $\sigma^{2}_{`s`}$) parameters.
  The filled circles represent the average $\mathrm{d}N/\mathrm{d}S$
  estimates while the empty squares represent the average $\omega$ estimates,
  across $20$ independent codon sequence alignments. The widths of the arrows
  correspond to twice the standard errors. The dashed lines highlight point
  of neutral selection effect.](FIG2.pdf){ width=38% }

Estimates of $\omega$ and $\mathrm{d}N/\mathrm{d}S$ summarised in Figure
\ref{testfig} strongly agree, except for the case of
$(\sigma^{2}_{`n`},\sigma^{2}_{`s`})=(0.10,0.02)$. The suppressed $\omega$
estimates, that is most pronounced for $\sigma^{2}_{`n`},\sigma^{2}_{`s`}>0$,
is likely a consequence of the well-documented conservative property of
homogeneous $\omega$ inference techniques (see for example, @nielsen1998).
Regardless, a correlation coefficient of approximately $0.9971$ was obtained
when the $\omega$ and $\mathrm{d}N/\mathrm{d}S$ averages were compared. The
standard errors ranged between $[0.0000,\;0.0077)$ and $(0.0004,\;0.0627)$
for the $\mathrm{d}N/\mathrm{d}S$ and $\omega$ estimates respectively. These
measures confirm that the outputs from `scoup` are accurate.

# Conclusions

We present `scoup`, a R package that allows for simulation of codon
sequences in a way that is capable of recapitulating the evolutionary
processes of biological systems more realistically than most existing
simulators. Our framework creatively incorporates the Ornstein-Uhlenbeck
process into the mutation-selection evolutionary model. This attribute
could potentially unlock exciting research avenues that will improve
existing knowledge about the complex interactions of different,
potentially interacting, molecular evolutionary processes. In another
unique contribution to the literature, the magnitude of the Darwinian
selection affect on the simulated sequences was controlled with the ratio
of the variances of selection coefficients. Given the summaries
in Figure \autoref{testfig}, we state the following hypothesis with respect
to natural selection inference from multi-population genetic sequences.
With $\omega$, it is difficult to fully distinguish between compensatory and
adaptive diversifying selection occurring on static and changing landscapes,
respectively. To establish this hypothesis, at least numerically,
`scoup` should be an invaluable resource.

# Code availability

`scoup` is published for free public use under the GPL-2 license. It
is available for download from the
[Bioconductor platform](doi.org/10.18129/B9.bioc.scoup), along with detailed
documentation and tutorial files. Additional help files may be accessed on
[Github](https://github.com/thsadiq/scoup).

# Whitepaper

The `scoup` whitepaper is available on the
[biorXiv](https://www.biorxiv.org/collection/evolutionary-biology) preprint
server.

# Acknowledgements

We thank Ben Murrell for suggesting modelling varying selection coefficients
with an OU process. Computations were performed using the
[`HPC1`](http://www.sun.ac.za/hpc) facility at Stellenbosch University, South
Africa.


# References

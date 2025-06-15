# flascs
Fitness Landscape Approach to Summarising Codon Sequences

## About
Curated functions that are designed to facilitate heuristic visual and quantitative summaries of protein sequences as part of mechanistic macro-evolutionary analyses. The aim is to demonstrate some interesting applications of the adaptive fitness landscape metaphor beyond the contexts of population genetic and quasispecies theories. The functions are targeted at genetic fitness landscapes. A primary objective is to extract some exploratory graphical and numerical details that will equip analysts with informed expectation of the evolutionary signatures present in input observed or simulated molecular sequences. Semi-theoretical and fully heuristic fitness landscape techniques are accessible, and they can be applied to sequences defined from any character space.

## Contents
- DESCRIPTION: file describing the package more elaborately.
- inst: a folder containing a temporary citation file
- man: a folder where the R markdown files are saved.
- NAMESPACE: what functions and/or classes are imported and/or exported.
- NEWS: Details of the edits, bug fixes and updates made to the package.
- R: a folder of R scripts that contains different (back-end) functions of the package.
- test: a directory of R scripts that contain some test functions. The contents are designed to facilitate bug detection.
- vignettes: a folder where the markdown and citation files for building the supporting tutorial document are saved.

## Installation
`flascs` can be installed from the Bioconductor platform with the following command. Interested users should visit the `devel` section on the Bioconductor website for access to the developmental version of the package that include the most recent updates and bug fixes.

``` r
if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
 BiocManager::install("flascs")
```

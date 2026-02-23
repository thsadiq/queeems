# queeems
Quantify the Extent of Evolutionary Evidence in Molecular Sequences

## About
Molecular protein sequences are only as useful, for evolutionary inferences, as the degree of signatures they retain. When very little time has passed since evolution from the most recent common ancestor, the information accrued in the sequences will most likely be insufficient for any valuable inference to be deducible. When the evolution time is too long, the sequences are said to have saturated. That is, the evolutionary signatures are most probably eroded for such sequences to yield reliable inferences. Functions designed to assist with quantifying the extent of saturation of protein sequences are presented. Compared to existing tools for similar purpose, the primary function herein add at least two new features to the literature. First, it is developed upon Bayesian principles. Second, it allows for site-wise assessment of saturation such that it useful for more than just phylogeny reconstruction inquests. For example, the package is amenable to heterogeneous natural selection pressure analyses.

## Contents
- [CONTRIBUTING](CONTRIBUTING.md): file where information about how to propose updates or seek assistance when struggling to execute any of the functions within the package.
- [DESCRIPTION](DESCRIPTION): file describing the package more elaborately.
- [inst](inst): a folder containing a temporary citation file and subfolder(s) where manuscripts of summarised description or of heavy applications of the package are saved.
- [LICENSE](LICENSE): a file that contains details of the GNU v3.0 that governs the free use and sharing of the package
- [man](man): a folder where the R markdown files that contain explanations of how to use each of the functions within the package are saved.
- [NAMESPACE](NAMESPACE): what functions and/or classes are imported and/or exported.
- [NEWS](NEWS): Details of the edits, bug fixes and updates made to the package.
- [R](R): a folder of R scripts that contains different (back-end) functions of the package.
- [README](README.md): the file that contains the package's preface information.
- [tests](tests): a directory of R scripts that contain some test functions. The contents are designed to facilitate internal bug detection.
- [vignettes](vignettes): a folder where the markdown and citation files for building the supporting tutorial document are saved.


## Installation
Upon completion `queeems` can be installed directly from GitHub with the following code.

``` r
if(!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
 remotes::install_github("thsadiq/queeems")
```

# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><       Reference Function & Necessary Classes and Generics.       ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         7th October 2025                         ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Background Functions and Data
aacids <- LETTERS[c(1,seq(3,9),seq(11,14),seq(16,20),22,23,25)]
aminoacid <- GENETIC_CODE[-which(GENETIC_CODE == "*")]
senseCodon <- sort( names(aminoacid))
nucleotides <- c("A", "C", "G", "T")

# ><>< # Define Classes
setClass("citing", representation(citeText="character"))
setClass("baseSummary", representation(baseType="character",
    extantSize="numeric", seqsNumerics="array"))
setClass("cncsframe", representation(nscensus="integer", wildseq="array",
    nstype="character", maxdiff="numeric", invariants="integer"))
setClass("siteindices", representation(eindices="numeric", noinfo="numeric",
    meth="character", alphaR="numeric"))
setClass("saturateBF", representation(bf.pv="numeric",
    nuLL="numeric", altLL="numeric", noInfo="numeric", nsites="numeric"))
setClass("geneindex", representation(eindex="numeric", noinfo="numeric",
    meth="character", alphaR="numeric", nbases="numeric", synonym="logical"),
    prototype(alphaR=as.numeric(NA), synonym=as.logical(NA)))

# ><>< # Specify Generics
setGeneric("bstringCodons",
    function(bstringset) standardGeneric("bstringCodons") )
setGeneric("codonDissimilarity",
    function(bstringset, diffUnit) standardGeneric("codonDissimilarity") )
setGeneric("nseqs", function(bsumm) standardGeneric("nseqs") )
setGeneric("datatype", function(bsumm) standardGeneric("datatype") )
setGeneric("basecensus", function(bsumm) standardGeneric("basecensus") )
setGeneric("nORs", function(cNS) standardGeneric("nORs") )
setGeneric("nsfreqs", function(cNS) standardGeneric("nsfreqs") )
setGeneric("nonvaries", function(cNS) standardGeneric("nonvaries") )
setGeneric("nucbalance", function(cNS) standardGeneric("nucbalance") )
setGeneric("wildtriplets", function(cNS) standardGeneric("wildtriplets") )
setGeneric("renyiA", function(sent) standardGeneric("renyiA") )
setGeneric("informula", function(sent) standardGeneric("informula") )
setGeneric("sitentropies", function(sent) standardGeneric("sitentropies") )
setGeneric("useSyn", function(gent) standardGeneric("useSyn") )
setGeneric("gentropy", function(gent) standardGeneric("gentropy") )
setGeneric("sitecount", function(gent) standardGeneric("sitecount") )
setGeneric("BFs", function(satube) standardGeneric("BFs") )
setGeneric("LogL0", function(satube) standardGeneric("LogL0") )
setGeneric("LogL1", function(satube) standardGeneric("LogL1") )

# ><>< # Create Citation Function
txt0 <- paste("Below is/are citation(s) relevant to the queeems package:")
txt1 <- paste("The corresponding BibTeX entry(ies) is/are as follows:")
tr1 <- paste0("Sadiq, H. (in progress). queeems: Quantify the Extent ",
    "of Evolutionary Evidence in Molecular Sequences. R package.")
tr2 <- paste0("Sadiq, H. (in progress). Saturation ",
    "Threshold of Codon Evolutionary Models. Preprint.")
br1 <- paste0("@Manual{,\n\t","title = {queeems: Quantify the Extent ",
    "of Evolutionary Evidence in Molecular Sequences},\n\t","author =",
    " {Hassan Sadiq},\n\tyear = {in progress},\n\t", "note = {R",
    " Package},\n\turl = {https://github.com/thsadiq/queeems},\n}")
br2 <- paste0("@Manual{,\n\t","title = {Saturation Threshold of Codon ",
    "Evolutionary Models},\n\t","author = {Hassan Sadiq},\n\t",
    "year = {in progress},\n\t", "note = {Preprint},\n}")
citevect <- c(tr1, tr2, br1, br2)
citeData <- matrix(citevect, nrow=length(citevect)%/%2, byrow=FALSE)

aboutQueeems <- function(cite=NULL){
citemsg <- paste0("Invalid citation index. There are ",
    nrow(citeData), " available queeems-related references.")
    input <- ifelse(is(cite,"NULL"), TRUE, cite %in% seq(1,nrow(citeData)))
    if(!input) stop( citemsg )
    if( is(cite,"NULL") ){
        citeID <- seq(1, nrow(citeData)) } else { citeID <- cite }
    requestedCite <- paste(c(txt0, citeData[citeID,1],
        txt1, citeData[citeID,2]), collapse="\n\n")
    output <- new("citing", citeText=requestedCite)
    return( output )
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

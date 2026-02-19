# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  Methods for `geneindex` Class.                  ><>< #
# ><>< ================================================================ ><>< #

setMethod("renyiA", "geneindex", function(sent) sent@alphaR )
setMethod("nonvaries", "geneindex", function(cNS) cNS@noinfo )
setMethod("informula", "geneindex", function(sent) sent@meth )
setMethod("useSyn", "geneindex", function(gent) gent@synonym )
setMethod("gentropy", "geneindex", function(gent) gent@eindex )
setMethod("sitecount", "geneindex", function(gent) gent@nbases )

gget <- function(gent){
    csyn <- useSyn(gent)
    hvalue <- gentropy( gent)
    colsize <- sitecount(gent)
    blanks <- length( nonvaries(gent))
    m0txt <- toupper( informula(gent))
    aword <- sprintf("%.2f", renyiA(gent))
    rtext <- paste0("Renyi (\U03B1 \U2248 ", aword,")")
    mword <- ifelse(m0txt=="SHANNON", "Shannon", rtext)
    sOut <- list(mword=mword, toutcol=colsize,
        void=blanks, hindex=hvalue, csyn=csyn)
    return(sOut)
}

setMethod("show", "geneindex", function(object) {
    ans0 <- gget(object)
    cat("\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    cat("\n\aSummarised Output for Entropy Analyses:")
    cat("\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    cat("\nEntropy Estimation Method:\t ", ans0$mword)
    if(!is(ans0$syn,"NA")) cat("\nSynonymous Counts:\t\t ", ans0$csyn)
    cat("\nNon-informative Sites:\t\t ", ans0$void, "/", ans0$toutcol)
    cat("\nProtein Entropy (H):\t\t ", sprintf("%.4f", ans0$hindex))
    cat("\n========================================================\n\n")
    }
)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

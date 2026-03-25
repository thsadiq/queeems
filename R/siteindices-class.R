# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                 Methods for `siteindices` Class.                 ><>< #
# ><>< ================================================================ ><>< #

setMethod("renyiA", "siteindices", function(sent) sent@alphaR )
setMethod("nonvaries", "siteindices", function(cNS) cNS@noinfo )
setMethod("informula", "siteindices", function(sent) sent@meth )
setMethod("sitentropies", "siteindices", function(sent) sent@eindices )

sgroup <- function(sent){
    avgH <- mean( sitentropies(sent))
    stale <- length( nonvaries(sent))
    colsize <- length( sitentropies(sent))
    serr <- sd( sitentropies(sent)) / colsize
    mOFe <- serr * sqrt(1 / 0.05)
    botci <- avgH - mOFe
    topci <- avgH + mOFe
    sOut <- list(void=stale, toutcol=colsize,
        expH=avgH, stderr=serr, lowci=botci, uppci=topci)
    return(sOut)
}

setMethod("show", "siteindices", function(object) {
    ans0 <- sgroup(object)
    m0sub <- toupper( informula(object))
    atext <- sprintf("%.2f", renyiA(object))
    redit <- paste("Renyi (\U03B1 \U2248 ", atext,")", sep="")
    m1sub <- ifelse(m0sub=="SHANNON", "Shannon", redit)
    cat("\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    cat("\n\aSummarised Output for Entropy Analyses:")
    cat("\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    cat("\nEntropy (H) Estimation Method:\t ", m1sub)
    cat("\nNon-informative Sites:\t\t ", ans0$void, "/", ans0$toutcol)
    cat("\nMean of H:\t\t\t ", sprintf("%.4f", ans0$expH))
    cat("\nStandard Error of H:\t\t ", sprintf("%.4f", ans0$stderr))
    cat("\nApprox. Lower 95% C.I. of H:\t ", sprintf("%.4f", ans0$lowci))
    cat("\nApprox. Upper 95% C.I. of H:\t ", sprintf("%.4f", ans0$uppci))
    cat("\n========================================================\n\n")
    }
)

setMethod("summary", "siteindices", function(object) {
    ans0 <- sgroup(object)
    relevant <- data.frame( Entropy=c(noinfo.sites=ans0$void,
        full.sites=ans0$toutcol, mean.H=ans0$expH, se.H=ans0$stderr,
        lower.ci=ans0$lowci, upper.ci=ans0$uppci))
    return(relevant)
    }
)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

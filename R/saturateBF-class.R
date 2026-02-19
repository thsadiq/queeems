# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  Methods for `saturateBF` Class                  ><>< #
# ><>< ================================================================ ><>< #

setMethod("BFs", "saturateBF", function(satube) satube@bf.pv )
setMethod("LogL0", "saturateBF", function(satube) satube@nuLL )
setMethod("nonvaries", "saturateBF", function(cNS) cNS@noInfo )
setMethod("LogL1", "saturateBF", function(satube) satube@altLL )
setMethod("sitecount", "saturateBF", function(gent) gent@nsites )

bfproc <- function(satube){
    h0limit <- 1
    h1limit <- 10
    ll0 <- LogL0(satube)
    ll1 <- LogL1(satube)
    gLL0 <- sum(ll0)
    gLL1 <- sum(ll1)
    homoBF <- gLL1 / gLL0
    nv <- sitecount(satube)
    bayesfacts <- BFs(satube)
    sngltn <- length( nonvaries(satube))
    satsites <- which(bayesfacts > h1limit)
    unclear <- (homoBF > h0limit) & (homoBF < h1limit)
    decide <- ifelse(homoBF > h1limit, "Not Saturated",
        ifelse(unclear, "Inconclusive", "Saturated"))
    sOut <- list(nLength=nv, fixed=sngltn, ll0=ll0,
        ll1=ll1, bf=bayesfacts, saturated=satsites,
        gH0=gLL0, gH1=gLL1, gBF=homoBF, gFin=decide)
    return(sOut)
}

showsout <- function(slist, app){
    detect <- length(slist$saturated)
    w0 <- "\nLogL(Saturated \U03B8\U2080 | Full Data):\t\t"
    w1 <- "\nLogL(Not Saturated \U03B8\U2081 | Full Data):\t"
    cat(paste0("\n",paste(rep(":",59),collapse="")))
    cat( paste0("\nSaturation Analyses Output Summary: ",app))
    cat(paste0("\n",paste(rep(":",59),collapse="")))
    cat("\nInvariant Sites:\t\t\t ", slist$fixed, "/", slist$nLength)
    cat("\nSaturated Sites:\t\t\t ", detect, "/", slist$nLength)
    cat(w0, sprintf("% .2f", slist$gH0))
    cat(w1, sprintf("% .2f", slist$gH1))
}

setMethod("show", "saturateBF", function(object) {
    ans0 <- bfproc(object)
    priorprint <- showsout(ans0, "Bayesian")
    cat("\nGlobal Bayes Factor:\t\t\t ", sprintf("%.2f", ans0$gBF))
    cat("\nGlobal Decision:\t\t\t ", ans0$gFin)
    cat(paste0("\n",paste(rep("=",59),collapse=""),"\n\n"))
    }
)

setMethod("summary", "saturateBF", function(object) {
    ans0 <- bfproc(object)
    relevant <- data.frame( Global.Summary=c(invary=ans0$fixed,
        nsites=ans0$nLength, satur.s=length(ans0$saturated),
        null.logL=ans0$gH0, alt.logL=ans0$gH1,
        fullBF=ans0$gBF, seq.state=ans0$gFin))
    return(relevant)
    }
)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

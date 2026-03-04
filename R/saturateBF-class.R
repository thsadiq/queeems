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
    h1limit <- 100
    ll0init <- LogL0(satube)
    ll1init <- LogL1(satube)
    gLL0 <- sum(ll0init)
    gLL1 <- sum(ll1init)
    homoBF <- ifelse(gLL1>0 & gLL0>0, gLL1/gLL0, (gLL1/gLL0)^(-1))
    nv <- sitecount(satube)
    bayesfacts <- BFs(satube)
    sngltn <- length( nonvaries(satube))
    satsites <- which(bayesfacts > h1limit)
    decide <- ifelse(homoBF < h1limit, "Not Saturated", "Saturated")
    sOut <- list(nLength=nv, fixed=sngltn, ll0=ll0init,
        ll1=ll1init, bf=bayesfacts, saturated=satsites,
        gH0=gLL0, gH1=gLL1, gBF=homoBF, gFin=decide)
    return(sOut)
}

showsout <- function(slist, app){
    detect <- length(slist$saturated)
    w0 <- "LogL(Not Saturated \U03B8\U2080 | Full Data):\t"
    w1 <- "LogL(Saturated \U03B8\U2081 | Full Data):\t\t"
    writeLines(paste0("\n\n",paste(rep(":",59),collapse="")))
    writeLines(paste0("Saturation Analyses Output Summary: ",app))
    writeLines(paste0("",paste(rep(":",59),collapse="")))
    writeLines(paste0("Invariant Sites:\t\t\t ",
        slist$fixed, "/", slist$nLength))
    writeLines(paste0("Saturated Sites:\t\t\t ",detect,"/",slist$nLength))
    writeLines(paste0(w0, sprintf("% .2f", slist$gH0)))
    writeLines(paste0(w1, sprintf("% .2f", slist$gH1)))
}

setMethod("show", "saturateBF", function(object) {
    ans0 <- bfproc(object)
    priorprint <- showsout(ans0, "Bayesian")
    cat("Overall Bayes Factor:\t\t\t", sprintf("%.2f", ans0$gBF))
    cat("\nOverall Decision:\t\t\t", ans0$gFin)
    cat(paste0("\n",paste(rep("=",59),collapse=""),"\n\n"))
    }
)

setMethod("summary", "saturateBF", function(object) {
    ans0 <- bfproc(object)
    relevant <- data.frame( Overall.Summary=c(no.info=ans0$fixed,
        nsites=ans0$nLength, satur.s=length(ans0$saturated),
        null.logL=round(ans0$gH0,4), alt.logL=round(ans0$gH1,4),
        fullBF=round(ans0$gBF,4), seq.state=ans0$gFin))
    return(relevant)
    }
)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

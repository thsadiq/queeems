# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  Methods for `saturatePV` Class                  ><>< #
# ><>< ================================================================ ><>< #

setMethod("PVs", "saturatePV", function(satupe) satupe@bf.pv )

pvproc <- function(satupe){
    ll0 <- LogL0(satupe)
    ll1 <- LogL1(satupe)
    gLL0 <- sum(ll0)
    gLL1 <- sum(ll1)
    nv <- sitecount(satupe)
    gcrit <- qchisq(0.95, 2)
    probvalues <- PVs(satupe)
    gstat <- 2 * (gLL1 - gLL0)
    sngltn <- length( nonvaries(satupe))
    satsites <- which(probvalues > 0.05)
    gpval <- round( 1 - pchisq(gstat,2), 4)
    decide <- ifelse(gpval < 0.05, "Not Saturated", "Saturated")
    sOut <- list(nLength=nv, fixed=sngltn, ll0=ll0, ll1=ll1,
        saturated=satsites, gH0=gLL0, gH1=gLL1, gstat=gstat,
        gcrit=gcrit, gpval=gpval, gFin=decide)
    return(sOut)
}

setMethod("show", "saturatePV", function(object) {
    ans0 <- pvproc(object)
    priorprint <- showsout(ans0, "LRT")
    cat("\nGlobal \U03C7\U00B2 test stat.:\t\t\t ",sprintf("%.2f",ans0$gstat))
    cat("\nGlobal \U03C7\U00B2 crit. value:\t\t\t ",sprintf("%.2f",ans0$gcrit))
    cat("\nOne-tail p-value:\t\t\t ", sprintf("%.4f", ans0$gpval))
    cat("\nGlobal Decision:\t\t\t ", ans0$gFin)
    cat(paste0("\n",paste(rep("=",59),collapse=""),"\n\n"))
    }
)

setMethod("summary", "saturatePV", function(object) {
    ans0 <- pvproc(object)
    relevant <- data.frame( Global.Summary=c(invary=ans0$fixed,
        nsites=ans0$nLength, satur.s=length(ans0$saturated),
        null.logL=ans0$gH0, alt.logL=ans0$gH1, test.stat=ans0$gstat,
        crit.value=ans0$gcrit, p.value=ans0$gpval, seq.state=ans0$gFin))
    return(relevant)
    }
)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

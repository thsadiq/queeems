# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  Methods for `saturateBF` Class                  ><>< #
# ><>< ================================================================ ><>< #

setMethod("BFs", "saturateBF", function(satube) satube@bf.pv )
setMethod("LogL0", "saturateBF", function(satube) satube@nuLL )
setMethod("nonvaries", "saturateBF", function(cNS) cNS@noInfo )
setMethod("LogL1", "saturateBF", function(satube) satube@altLL )
setMethod("mainPi", "saturateBF", function(satube) satube@sspi )
setMethod("sitecount", "saturateBF", function(gent) gent@nsites )
setMethod("BFthreshold", "saturateBF", function(satube) satube@bfcutoff )

bfproc <- function(satube){
    h1limit <- BFthreshold(satube)
    ll0init <- LogL0(satube)
    ll1init <- LogL1(satube)
    nv <- sitecount(satube)
    sspi <- mainPi(satube)
    sspiedit <- (ceiling(1 / sspi) * 2) + 2
    sspiN <- max(102, sspiedit)
    bayesfacts <- BFs(satube)
    satlogic <- bayesfacts > h1limit
    sngltn <- length( nonvaries(satube))
    satPI <- seq(0, 1, length=sspiN)[seq(2,sspiN-1)]
    infcounts <- c(sum(satlogic), sum(!satlogic))
    globaLikely <- vapply(seq(1,sspiN-2), function(i0){
        mnomLogl(infcounts, c(satPI[i0],1-satPI[i0]))},
        FUN.VALUE=vector("numeric",1))
    globalPosterior <- softmax(globaLikely, TRUE)
    gALT <- globalPosterior[ which(satPI > mainPi(satube))]
    gNULL <- globalPosterior[ which(satPI <= mainPi(satube))]
    homoBF <- max(gALT) / max(gNULL)
    decide <- ifelse(homoBF <= h1limit, "Not Saturated", "Saturated")
    sOut <- list(nLength=nv, fixed=sngltn, ll0=ll0init, ll1=ll1init,
        bf=bayesfacts, saturated=which(satlogic), gBF=homoBF, gFin=decide)
    return(sOut)
}

showsout <- function(slist, app){
    detect <- length(slist$saturated)
    writeLines(paste("\n\n",paste(rep(":",59),collapse=""),sep=""))
    writeLines(paste("Saturation Analyses Output Summary: ",app,sep=""))
    writeLines(paste("",paste(rep(":",59),collapse=""),sep=""))
    writeLines(paste("Invariant Sites:\t\t\t ",
        slist$fixed, "/", slist$nLength,sep=""))
    writeLines(paste("Saturated Sites:\t\t\t ",
        detect,"/",slist$nLength,sep=""))
}

setMethod("show", "saturateBF", function(object) {
    ans0 <- bfproc(object)
    priorprint <- showsout(ans0, "Bayesian")
    cat("Overall Bayes Factor:\t\t\t", sprintf("%.2f", ans0$gBF))
    cat("\nOverall Decision:\t\t\t", ans0$gFin)
    cat(paste("\n",paste(rep("=",59),collapse=""),"\n\n",sep=""))
    }
)

setMethod("summary", "saturateBF", function(object) {
    ans0 <- bfproc(object)
    relevant <- data.frame( Overall.Summary=c(no.info=ans0$fixed,
        nsites=ans0$nLength, satur.s=length(ans0$saturated),
        fullBF=round(ans0$gBF,4), seq.state=ans0$gFin))
    return(relevant)
    }
)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

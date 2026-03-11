# ><>< ================================================================= ><>< #
# ><><                     Protein Sequence Analyses                     ><>< #
# ><><                             ~~~~~~~~~                             ><>< #
# ><><        Changing selection pressure analyses for vignettes.        ><>< #
# ><>< ================================================================= ><>< #

# ><>< # Create useful summary function
infoests <- function(xdata){
    serr <- sd(xdata) / sqrt( length(xdata))
    avg <- mean(xdata)
    low <- avg  -  serr
    upp <- avg  +  serr
    outs <- c(low=low, avg=avg, upp=upp)
    return( outs )
}

# ><>< # Function that makes plotting easier
tfunc <- function(i, pty, clr, cx, mat, agl=90){
    arrows(i, mat["low",i], i, mat["upp",i], .1, agl, 3, colors()[clr], lwd=3)
    points(i,mat["avg",i],pch=pty,lwd=2,col=colors()[clr],cex=cx,bg="white")
    return(0)
}

# ><>< # Load package into current R session
library("queeems")

# ><>< # Indicate number of simulated data replicates
nreps <- 5

# ><>< # Set preferred Bayes factor threshold
bflimit <- 3

# ><>< # Set preferred saturated sites tolerance
ssprop <- 0.01

# ><>< # Retrieve path to simulated tutorial data folder
datadir <- system.file("extdata/example/varNS/", package="queeems")

# Non-synonymous variance value
nsvary <- seq(1, 5) / 10
ntags <- sprintf("%.0f", nsvary*10)

# ><>< # Create data frames to save interesting outputs
cnames <- paste0("ns", ntags)
rnames <- c("low","avg","upp")
homoBF <- matrix(NA, nreps, length(ntags), dimnames=list(NULL,cnames))
sbayes <- matrix(NA, nreps, length(ntags), dimnames=list(NULL,cnames))

for(a1 in seq(1,length(ntags))){
    # ><>< # Read merged molecular sequences data as text
    seqsname <- file.path(datadir, paste0("gdata",ntags[a1],".txt"))
    seqstext <- readLines(seqsname)

    # ><>< # Extract information about number of simulated replicates
    sbreaks <- which(seqstext == "  64  300")

    # ><>< # Analyse data replicates captured from merged data file
    for(a0 in seq(1,nreps)){

        # ><>< # Extract and temporarily save a sequence replicate
        init <- sbreaks[a0] + 1
        halt <- sbreaks[a0] + 128
        seqrep <- seqstext[seq(init,halt)]
        tempfile <- file.path(tempdir(), "iseqs.fasta")
        write.table(seqrep, file=tempfile, append=FALSE,
            quote=FALSE, row.names=FALSE, col.names=FALSE)

        # ><>< # Execute saturation analysis
        satest <- seqSaturation(tempfile, "A", 6, 0.4, 7, bflimit, ssprop)

        # ><>< # Extract Bayes factors
        homoBF[a0,a1] <- as.numeric( summary(satest)["fullBF",])
        sbayes[a0,a1] <- sum(BFs(satest) > bflimit)
    }
    # ><>< # Provide analysis progress update
    message("Completed analyses for \U03c3\U00b2 = 0.", ntags[a1])
}
writeLines("\n\n# ><>< # Overall log(Bayes factor) estimates")
print( log(homoBF))

# ><>< # Summarise overall Bayes factor and dN/dS values
bfsumm <- apply(log(homoBF),2,infoests)
dndsPath <- file.path(datadir, "dnds.csv")
dndses <- read.table(dndsPath, sep=";", header=TRUE)

# ><>< # Facilitate super-imposition of dN/dS and BFs
transplot <- function(x, xrng, newrng){
    xmin <- xrng[1]
    ymin <- newrng[1]
    xdiff <- diff(xrng)
    ydiff <- diff(newrng)
    newx <- (((x - xmin) / xdiff) * ydiff) + ymin
    return(newx)
}
newaxs <- seq(4,9)
oldaxs <- transplot(newaxs, c(5,9), c(-.4,2.3))
newDnDs <- transplot(dndses, c(5,9), c(-.4,2.3))
dndsumm <- apply(newDnDs, 2, infoests)

# ><>< # Add summaries of the BF values to plot
par(mar=c(4.4,4.2,0.2,3.0))
plot(0, 0, col="white", bty="n", xlab="", ylab="",
    xlim=c(1.0,5.0), ylim=c(-.4,2.3), xaxt="n", yaxt="n")
phold <- vapply(seq(1,5), tfunc, 15, 50, 2.5, bfsumm, FUN.VALUE=0)
axis(1, seq(1,5), nsvary, tck=-.01, lwd=2, padj=-.5, cex.axis=1.5)
axis(2, las=1, tck=-0.01, lwd=2, cex.axis=1.5, hadj=0.7)
mtext(expression(sigma[n]^2), 1, 3.4, cex=2)
mtext("Log(Bayes factor)", 2, 2.5, cex=2)
text(1.2, 2.2, "B", cex=2.5, font=2)

# ><>< # Include summaries of the dN/dS values
phold <- vapply(seq(1,5), tfunc, 23, 17, 2, dndsumm, 60, FUN.VALUE=0)
axis(4, oldaxs, newaxs, las=1, tck=-0.01, lwd=2, cex.axis=1.5, hadj=0.7)
mtext("dN/dS", 4, 2.0, cex=2)

# ><>< # Add legend to plot
text(c(2,2), c(2.2,2.0), c("BF","dN/dS"), cex=1.8, pos=4)
arrows(1.7, 2.21, 2.0, 2.21, .1, 90, 3, colors()[50], 1, 3)
arrows(1.7, 2.0, 2.0, 2.0, .1, 60, 3, colors()[17], 1, 3)
points(c(1.85,1.85), c(2.21,2.0), cex=2, bg="white",
    lwd=3, pch=c(15,23), col=colors()[c(50,17)])
box(lwd=2)

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #

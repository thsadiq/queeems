# ><>< ================================================================= ><>< #
# ><><                     Protein Sequence Analyses                     ><>< #
# ><><                             ~~~~~~~~~                             ><>< #
# ><><            Changing stem length analyses for vignettes            ><>< #
# ><>< ================================================================= ><>< #

# ><>< # Create useful summary function
infoests <- function(xdata){
    serr <- sd(xdata) / sqrt( length(xdata))
    avg <- mean( xdata )
    low <- avg  -  serr
    upp <- avg  +  serr
    outs <- c(low=low, avg=avg, upp=upp)
    return( outs )
}

# ><>< # Load package into current R session
library("queeems")

# ><>< # Retrieve path to simulated tutorial data folder
datadir <- system.file("extdata/example/stemL/", package="queeems")

# ><>< # Indicate number of simulated data replicates
nreps <- 5

# ><>< # Set preferred Bayes factor threshold
bflimit <- 3

# ><>< # Set preferred saturated sites tolerance
ssprop <- 0.01

# Branch length
blent <- seq(0.03, 0.15, 0.03)
btags <- sprintf("%02.0f", blent*100)

# ><>< # Create data frames to save interesting outputs
rnames <- c("low","avg","upp")
cnames <- paste("bl", btags, sep="")
vnames <- paste("rep.", seq(1,nreps), sep="")
homoBF <- matrix(NA, nreps, length(btags), dimnames=list(vnames,cnames))
sbayes <- matrix(NA, nreps, length(btags), dimnames=list(vnames,cnames))
smrybx <- matrix(NA, 3, length(btags), dimnames=list(rnames,cnames))
cdnTmp <- aasTmp <- vector("numeric", 100*nreps)
nucTmp <- vector("numeric", 300*nreps)
nucent <- codent <- aasent <- smrybx

for(a1 in seq(1,length(btags))){
    # ><>< # Read merged molecular sequences data as text
    seqsname <- file.path(datadir, paste("seqs",btags[a1],".txt",sep=""))
    seqstext <- readLines(seqsname)

    # ><>< # Extract information about number of simulated replicates
    sbreaks <- which(seqstext == "  64  300")

    # ><>< # Analyse data replicates captured from merged data file
    for(a0 in seq(1,nreps)){
        nucentID <- seq(1+(a0-1)*300, a0*300)
        entropyID <- seq(1+(a0-1)*100, a0*100)

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

        # ><>< # Execute codon entropy analysis
        tempentropy <- molentropy(tempfile, "codon", "shannon")
        cdnTmp[entropyID] <- sitentropies(tempentropy)

        # ><>< # Process nucleotide entropy analysis
        tempentropy <- molentropy(tempfile, "dna", "shannon")
        nucTmp[nucentID] <- sitentropies(tempentropy)

        # ><>< # Transform nucleotide to amino acid sequences
        dnaformat <- Biostrings::readDNAStringSet(tempfile)
        aaformat <- Biostrings::translate(dnaformat)
        Biostrings::writeXStringSet(aaformat, tempfile)

        # ><>< # Implement amino acid entropy analysis
        tempentropy <- molentropy(tempfile, "aa", "shannon")
        aasTmp[entropyID] <- sitentropies(tempentropy)
    }
    # ><>< # Summarise and save entropy values
    nucent[rnames,a1] <- infoests( nucTmp )[rnames]
    codent[rnames,a1] <- infoests( cdnTmp )[rnames]
    aasent[rnames,a1] <- infoests( aasTmp )[rnames]

    # ><>< # Provide analysis progress update
    message("Completed analyses of stem length 0.", btags[a1])
}
writeLines("\n\n# ><>< # Site-wise Bayes factors peak")
print( head( BFs(satest)))

writeLines("\n\n# ><>< # Saturated site counts")
print(sbayes)

writeLines("\n\n# ><>< # Overall Bayes factor estimates")
print(homoBF)

# ><>< # Function that makes plotting easier
tfunc <- function(i, pty, clr, cx, mat, agl=90){
    arrows(i, mat["low",i], i, mat["upp",i], .1, agl, 3, colors()[clr], lwd=3)
    points(i,mat["avg",i],pch=pty,lwd=2,col=colors()[clr],cex=cx,bg="white")
    return(0)
}

# ><>< # Summarise overall Bayes factor values
bfsumm <- apply(homoBF,2,infoests)

par(mar=c(4.2,4.2,0.2,0.2))
plot(0, 0, col="white", bty="n", xlab="", ylab="",
    xlim=c(1.0,5.0), ylim=c(0.0,7.0), xaxt="n", yaxt="n")
phold <- vapply(seq(1,5), tfunc, 15, 50, 2.5, bfsumm, FUN.VALUE=0)
axis(1, seq(1,5), seq(3,15,3)/100*7, tck=-.01, lwd=2, padj=-.2, cex.axis=1.5)
axis(2, las=1, tck=-0.01, lwd=2, cex.axis=1.5, hadj=0.2)
text(1.3, 6.5, "A", cex=2.5, font=2)
mtext("Bayes factor", 2, 2.2, cex=2)
mtext("Tree Length", 1, 2.4, cex=2)
box(lwd=2)

# ><>< # Scale entropy estimates by the maximum value
eScaledNuc <- nucent / (-log(1/4))
eScaledCodon <- codent / (-log(1/61))
eScaledAmino <- aasent / (-log(1/20))

# ><>< # Plot distributions of entropy estimates
par(mar=c(4.2,4.2,0.2,0.2))
bases <- c("Nucleotide","Codon","Amino acid")
plot(0, 0, col="white", bty="n", xlab="", ylab="",
    xlim=c(1.0,5.0), ylim=c(.19,0.57), xaxt="n", yaxt="n")
phold <- vapply(seq(1,5), tfunc, 0, 17, 1.5, eScaledNuc, FUN.VALUE=0)
phold <- vapply(seq(1,5), tfunc, 1, 60, 1.5, eScaledCodon, FUN.VALUE=0)
phold <- vapply(seq(1,5), tfunc, 8, 74, 1.5, eScaledAmino, FUN.VALUE=0)
axis(1, seq(1,5), seq(3,15,3)/100*7, tck=-.01, lwd=2, padj=-.5, cex.axis=1.5)
text(rep(3.7,3), c(.27,.24,.21)-.002, bases, cex=1.5, pos=4)
axis(2, las=1, tck=-0.01, lwd=2, cex.axis=1.5, hadj=0.7)
points(rep(3.6,3), c(.27,.24,.21), cex=2.5, lwd=2,
    pch=c(0,1,8), col=colors()[c(17,60,74)])
text(3.8, 0.30, "Base", cex=2, font=2)
mtext("Scaled Entropy", 2, 2.4, cex=2)
text(1.3, 0.54, "A.1", cex=2.5, font=2)
mtext("Tree Length", 1, 2.2, cex=2)
box(lwd=2)

# ><>< # Site-wise Bayes factors from replicate 5 of length 1.05
par(mar=c(4.2,4.2,0.2,0.2))
plot(0, 0, col="white", bty="n", xlab="", ylab="",
    xlim=c(.6,4.4), ylim=c(.03,0.85), xaxt="n", yaxt="n")
mtext(c("Site-Wise Bayes Factors","Density"), c(1,2), 2.4, cex=2)
axis(2, las=1, tck=-0.01, lwd=2, cex.axis=1.5, hadj=0.7)
hist(BFs(satest), freq=FALSE, add=TRUE, xaxt="n",
    yaxt="n", main="", col=colors()[17], lwd=2)
axis(1, tck=-.01, lwd=2, padj=-.5, cex.axis=1.5)
text(4.1, 0.81, "A.2", cex=2, font=2)
box(lwd=2)


# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #

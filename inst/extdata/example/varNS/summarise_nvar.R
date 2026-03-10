# ><>< ================================================================= ><>< #
# ><><                     Protein Sequence Analyses                     ><>< #
# ><><                             ~~~~~~~~~                             ><>< #
# ><><            Changing stem length analyses for vignettes            ><>< #
# ><>< ================================================================= ><>< #

# ><>< # Create useful summary function
infoests <- function(xdata){
    serr <- sd(xdata) / sqrt( length(xdata))
    low <- avg  -  serr
    upp <- avg  +  serr
    avg <- mean(xdata)
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
cnames <- paste0("bl", btags)
rnames <- c("low","avg","upp")
homoBF <- matrix(NA, nreps, length(btags), dimnames=list(NULL,cnames))
sbayes <- matrix(NA, nreps, length(btags), dimnames=list(NULL,cnames))
smrybx <- matrix(NA, 3, length(btags), dimnames=list(rnames,cnames))
cdnTmp <- aasTmp <- vector("numeric", 100*nreps)
nucTmp <- vector("numeric", 300*nreps)
nucent <- codent <- aasent <- smrybx

for(a1 in seq(1,length(btags))){
    # ><>< # Read merged molecular sequences data as text
    seqsname <- file.path(datadir, paste0("seqs",btags[a1],".txt"))
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
    message("Completed analyses of stem length ", btags[a1])
}
writeLines("\n\n# ><>< # Saturated site counts")
print(sbayes)

writeLines("\n\n# ><>< # Overall Bayes factor estimates")
print(homoBF)

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #

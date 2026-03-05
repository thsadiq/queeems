# ><>< ================================================================= ><>< #
# ><><                     Protein Sequence Analyses                     ><>< #
# ><><                             ~~~~~~~~~                             ><>< #
# ><><  R code used to generate the results presented in the vignettes.  ><>< #
# ><>< ================================================================= ><>< #

rm(list=ls())
codePath <- "/Users/buyer/Dropbox/2024/ou-paper/bioc/queeems"
devtools::load_all( file.path(codePath, "R"))

# ><>< ==================== ><>< #

# ><>< # Load package into current R session
# library("queeems")

# ><>< # Retrieve path to simulated tutorial data folder
datadir <- system.file("extdata/example/", package="queeems")

# ><>< # Read merged molecular sequences data as text
seqsname <- file.path(datadir, "tut_seqs.txt")
seqstext <- readLines(seqsname)

# ><>< # Extract information about number of simulated replicates
sbreaks <- which(seqstext == "  64  750")
nreps <- length(sbreaks)

# ><>< # Create data frames to save interesting outputs
homoBF <- vector("numeric", nreps)
cnames <- paste0("r", seq(1,nreps))
nucent <- matrix(NA, 750, nreps, dimnames=list(NULL,cnames))
codent <- matrix(NA, 250, nreps, dimnames=list(NULL,cnames))
aasent <- matrix(NA, 250, nreps, dimnames=list(NULL,cnames))
sbayes <- matrix(NA, 250, nreps, dimnames=list(NULL,cnames))

# ><>< # Storage for saturation-free data to allow subsequent analyses
nosatpath <- file.path(tempdir(), "satfree.txt")

# ><>< # Adhoc function to extract saturation-free data
printsite <- function(b0, sqarray, ipath=nosatpath){
    ntag <- rownames(sqarray)[b0]
    sqline <- paste(sqarray[b0,], collapse="")
    fulltxt <- paste0(">", ntag, "\n", sqline)
    write.table(fulltxt, ipath, TRUE, FALSE, col.name=FALSE, row.names=FALSE)
    return(0)
}

# ><>< # Analyse data replicates captured from merged data file
for(a0 in seq(1,nreps)){

    # ><>< # Extract and temporarily save a sequence replicate
    init <- sbreaks[a0] + 1
    halt <- sbreaks[a0] + 128    # >< # Each replicate is from 64 leaves
    seqrep <- seqstext[seq(init,halt)]
    tempfile <- file.path(tempdir(), "iseqs.fasta")
    write.table(seqrep, file=tempfile, append=FALSE,
        quote=FALSE, row.names=FALSE, col.names=FALSE)

    # ><>< # Execute saturation analysis
    satest <- seqSaturation(tempfile, "A", 6, 0.4, 7)
    homoBF[a0] <- as.numeric( summary(satest)["fullBF",])
    sbayes[,a0] <- BFs(satest)

    # ><>< # Extract and save saturation-free sequence locations
    nullsites <- which(BFs(satest) < 100)
    sieved <- seqfilter(tempfile, nullsites, "codon")
    imatrix <- as.matrix(sieved)
    itext <- paste0("  ", nrow(imatrix), "  ", ncol(imatrix))
    write.table(itext,nosatpath,a0!=1,FALSE,row.names=FALSE,col.names=FALSE)
    pholder <- vapply(seq(1,nrow(imatrix)), printsite, imatrix, FUN.VALUE=0)
    write.table("\n",nosatpath,TRUE,FALSE,row.names=FALSE,col.names=FALSE)

    # ><>< # Execute codon entropy analysis
    entcdn <- molentropy(tempfile, "codon", "shannon")
    codent[,a0] <- sitentropies(entcdn)

    # ><>< # Process nucleotide entropy analysis
    entnuc <- molentropy(tempfile, "dna", "shannon")
    nucent[,a0] <- sitentropies(entnuc)

    # ><>< # Transform nucleotide to amino acid sequences
    dnaformat <- Biostrings::readDNAStringSet(tempfile)
    aaformat <- Biostrings::translate(dnaformat)
    Biostrings::writeXStringSet(aaformat, tempfile)

    # ><>< # Implement amino acid entropy analysis
    entamino <- molentropy(tempfile, "aa", "shannon")
    aasent[,a0] <- sitentropies(entamino)

    # ><>< # Provide analysis progress update
    message("Completed analyses of data replicate ", a0)
}

# ><>< # Read dN/dS values from file
dndsfile <- file.path(datadir, "tut_dnds.csv")
dndsvalues <- read.table(dndsfile, header=TRUE)

# ><>< # Read MLE of non-synonymous/synonymous substitution from file
mlefile <- file.path(datadir, "tut_mles.txt")
mlevalues <- read.table(mlefile, header=TRUE)

# ><>< # Temporary storage of the saturation free sequences.
# ><>< # Data was retrieved from where it is saved and analysed with PAML
message("\nSequences that are rid of saturated sites",
    " may be retrieved from:\n\t", nosatpath, "\n")

# ><>< # Generate distributions of entropy estimates
par(fig=c(0, 0.36, 0, 1), mar=c(4.2,4.2,0.2,0))
    plot(0, 0, col="white", bty="n", xlab="", ylab="",
        xlim=c(0,1.4), ylim=c(0.04,1.2), xaxt="n", yaxt="n")
    axis(1, tck=-0.015, lwd=2, padj=-.3, cex.axis=1.5)
    axis(2, las=1, tck=-0.015, lwd=2, cex.axis=1.5, hadj=0.7)
    a1 <- apply(nucent, 2, function(i) lines(density(i), lwd=.5, col="grey60"))
    lines(rep(-log(1/4),2), c(-0.1,1.3), lwd=2, lty=2, col="tomato")
    mtext("Nucleotide Entropy", 1, 2.7, cex=1.5)
    mtext("Density", 2, 2.4, cex=1.5)
    box(lwd=2)
par(fig=c(0.36, 0.68, 0, 1), mar=c(4.2,0,0.2,0), new=TRUE)
    plot(0, 0, col="white", bty="n", xlab="", ylab="",
        xlim=c(0,4.2), ylim=c(0.04,1.2), xaxt="n", yaxt="n")
    axis(1, tck=-0.015, lwd=2, padj=-.3, cex.axis=1.5)
    a1 <- apply(codent, 2, function(i) lines(density(i), lwd=.5, col="grey60"))
    lines(rep(-log(1/61),2), c(-0.1,1.3), lwd=2, lty=2, col="tomato")
    mtext("Codon Entropy", 1, 2.7, cex=1.5)
    box(lwd=2)
par(fig=c(0.68, 1, 0, 1), mar=c(4.2,0,0.2,0.2), new=TRUE)
    plot(0, 0, col="white", bty="n", xlab="", ylab="",
        xlim=c(0,3.1), ylim=c(0.04,1.2), xaxt="n", yaxt="n")
    axis(1, tck=-0.015, lwd=2, padj=-.3, cex.axis=1.5)
    a1 <- apply(aasent, 2, function(i) lines(density(i), lwd=.5, col="grey60"))
    lines(rep(-log(1/20),2), c(-0.1,1.3), lwd=2, lty=2, col="tomato")
    mtext("Amino Acid Entropy", 1, 2.7, cex=1.5)
    box(lwd=2)






# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #

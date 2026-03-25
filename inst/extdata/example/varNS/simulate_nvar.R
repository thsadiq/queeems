# ><>< ================================================================= ><>< #
# ><><                    Protein Sequence Simulation                    ><>< #
# ><><                             ~~~~~~~~~                             ><>< #
# ><><  R code used to simulate the protein sequences discussed in the   ><>< #
# ><><                analyses presented in the vignettes                ><>< #
# ><>< ================================================================= ><>< #

# ><>< # ======
# Combine sequences from >1 simulation executions to facilitate PAML analyses.
# - scoupObj: output of 'alignsim' from each simulation execution.
# - simID: index of the current execution i.e. order of alignment.
# - seqLoc: path to file where the merged data shoud be stored.
# ><>< # ======
mergeseq <- function(scoupObj, simID, seqLoc){
    newData <- cseq(scoupObj)
    nleaf <- nrow(newData)
    ndna <- ncol( seqs(scoupObj)) * 3
    apend <- ifelse(simID == 1, FALSE, TRUE)
    preSEQ <- paste("  ", nleaf, "  ", ndna, sep="")
    write.table(preSEQ, seqLoc, apend, FALSE,
            row.names=FALSE, col.names=FALSE)
    wrt <- vapply(seq(1,nleaf), function(w){
        sq0 <- paste(">S", sprintf("%03.0f",w), "\n", newData[w,], sep="")
        write.table(sq0,seqLoc,TRUE,FALSE,row.names=FALSE,col.names=FALSE)
        return(0)
        }, FUN.VALUE=0)
    write.table("\n", seqLoc, TRUE, FALSE, row.names=FALSE, col.names=FALSE)
    dndsVALUE <- mean( dNdS(scoupObj) )
    return(dndsVALUE)
}
# ><>< # ======

# Make simulation package accessible in R session
library(scoup)

# Number of tree leaves
leaves <- 64

# Number of codon sites
sitesize <- 100

# Number of data replications for each tree length
sims <- 5

# Branch length
blent <- 0.10

# Synonymous variance value
synvary <- 0.01

# Effective population size
effpop <- 1000

# OU reversion parameter (Theta) value
eThta <- 0.0

# OU asymptotic variance value
eVary <- 0.0

# Create appropriate simulation function ("ou") object
adaptStat <- ouInput(c(eVar=eVary,Theta=eThta))

# Sequence alignment size information
seqStat <- seqDetails(c(nsite=sitesize, ntaxa=leaves, blength=blent))

# Specify data storage paths
pryPath <- tempdir()

# Non-synonymous variance value
nsvary <- seq(1, 5) / 10
ntags <- sprintf("%.0f", nsvary*10)

# Create storage for dN/dS estimates
dndsPath <- file.path(pryPath, "dnds.csv")
dndsarray <- matrix(NA, ncol=length(nsvary),
        nrow=sims, dimnames=list(NULL,paste("ns0",ntags,sep="")))

for(h in seq(1,length(nsvary))){
    # OU landscape shift parameters
    newvnvs <- ifelse(synvary == 0, 0, nsvary[h]/synvary)
    hbrunoStat <- hbInput(c(vNvS=newvnvs, nsynVar=nsvary[h], Ne=effpop))

    # Create sequence file
    seqPath <- file.path(pryPath, paste("gdata",ntags[h],".txt",sep=""))

    # Iterate over the specified number of replicates
    for(i in seq(1,sims)){ 

        # Execute simulation
        simData <- alignsim(adaptStat, seqStat, hbrunoStat, NA)
        dndsavg <- mergeseq(simData, i, seqPath)
        dndsarray[i,h] <- dndsavg
    }
    message("Generated sequences for \U03c3\U00b2 = 0.", ntags[h])
}
write.table(dndsarray, dndsPath, FALSE, FALSE, ";", row.names=FALSE)
message("\nscoup simulation completed and outputs",
    " are saved in:\n\t", pryPath, "\n")

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #

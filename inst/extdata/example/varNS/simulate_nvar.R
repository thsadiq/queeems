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
    preSEQ <- paste0("  ", nleaf, "  ", ndna)
    write.table(preSEQ, seqLoc, apend, FALSE,
            row.names=FALSE, col.names=FALSE)
    wrt <- vapply(seq(1,nleaf), function(w){
        sq0 <- paste0(">S", sprintf("%03.0f",w), "\n", newData[w,])
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
blent <- seq(0.03, 0.15, 0.03)
btags <- sprintf("%02.0f", blent*100)

# Non-synonymous variance value
nsynvary <- 0.25

# Synonymous variance value
synvary <- 0.01

# Effective population size
effpop <- 1000

# OU landscape shift parameters
newvnvs <- ifelse(synvary == 0, 0, nsynvary/synvary)
hbrunoStat <- hbInput(c(vNvS=newvnvs, nsynVar=nsynvary, Ne=effpop))

# OU reversion parameter (Theta) value
eThta <- 0.0

# OU asymptotic variance value
eVary <- 0.0

# Create appropriate simulation function ("ou") object
adaptStat <- ouInput(c(eVar=eVary,Theta=eThta))

# Specify data storage paths
pryPath <- "stemL" # tempdir()

# Create storage for dN/dS estimates
dndsPath <- file.path(pryPath, "dnds.csv")
dndsarray <- matrix(NA,ncol=length(blent),nrow=sims,dimnames=list(NULL,btags))

for(h in seq(1,length(blent))){
    # Save the evolutionary tree used for the simulation
    t3newick <- biTree(leaves, blent[h])
    t3Path <- file.path(pryPath, paste0("tree", btags[h],".txt"))
    write.table(t3newick,t3Path,FALSE,FALSE,row.names=FALSE,col.names=FALSE)

    # Create sequence file
    seqPath <- file.path(pryPath, paste0("gdata", btags[h],".txt"))

    # Sequence alignment size information
    seqStat <- seqDetails(c(nsite=sitesize, ntaxa=leaves, blength=blent[h]))

    # Iterate over the specified number of replicates
    for(i in seq(1,sims)){ 

        # Execute simulation
        simData <- alignsim(adaptStat, seqStat, hbrunoStat, NA)
        dndsavg <- mergeseq(simData, i, seqPath)
        dndsarray[i,h] <- dndsavg
    }
    message("Generated sequences for stem length ", sprintf("%.02f",blent[h]))
}
write.table(dndsarray, dndsPath, FALSE, FALSE, ";", row.names=FALSE)
message("\nscoup simulation completed and outputs",
    " are saved in:\n\t", pryPath, "\n")

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #

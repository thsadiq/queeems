# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                          FUBAR Weights.                          ><>< #
# ><>< ================================================================ ><>< #

pispecs <- function(nweights, wless1, weightmax){
    wrd1 <- "\nNumber of weight categories (`nweights`) must be \U2265 4."
    wrd2 <- paste("\nInappropriate `wless1` input:\nMust be such that",
        " at least two weight values are > 0 and \U2264 1.", sep="")
    wrd3 <- paste("\nIncorrect weight upper bound input:",
                "\nRequire 'weightmax' > 1.", sep="")
    if(nweights <= 3) stop(wrd1, call.=FALSE)
    if(weightmax <= 1) stop(wrd3, call.=FALSE)
    if( floor(wless1 * nweights) < 2) stop(wrd2, call.=FALSE)
    return( 0 )
}

fubarweights <- function(nweights, wless1, weightmax, addzero){
    preamble <- pispecs(nweights, wless1, weightmax)
    fkappa <- floor(wless1 * nweights)
    fdeltatop <- (weightmax - 1)^(1/3)
    c1init <- ifelse(addzero, 0, 1)
    fdeltabot <- nweights - fkappa - (1 - c1init)
    fdelta <- fdeltatop / fdeltabot
    class1 <- seq(c1init, fkappa) / fkappa
    class2 <- 1 + ((seq(1, fdeltabot) * fdelta)^3)
    weightvector <- c(class1, class2)
    return(weightvector)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

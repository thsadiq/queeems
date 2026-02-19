# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         SoftMax Function                         ><>< #
# ><>< ================================================================ ><>< #

wcheck <- function(weights){
    shout0 <- "\nInvalid weights input. Check that input not all NA or zeros."
    shout1 <- "\n`NA`s detected in input weights replaced with zero(s)."
    rule1 <- all( is.na(weights)) | all(weights==0) | length(weights)==0
    if(rule1) stop(shout0, call.=FALSE)
    someNA <- is.na( weights)
    if(sum(someNA) > 0) warning(shout1, call.=FALSE)
    weights[someNA] <- 0
    return(weights)
}

softmax <- function(weights, maxadj){
    neWeights <- wcheck(weights)
    smax <- ifelse(maxadj, max(neWeights), 0)
    wShifted <- neWeights - smax
    expweights <- exp( wShifted )
    newFreq <- expweights / sum(expweights)
    return(newFreq)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

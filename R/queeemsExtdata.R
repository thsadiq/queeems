# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><             Extract the Paths to Example FASTA Data.             ><>< #
# ><>< ================================================================ ><>< #

queeemsExtdata <- function(filename=NULL){
    if(is(filename, "NULL")){
        filePath <- dir(system.file("extdata", package="queeems"))
    }else{
        filePath <- system.file("extdata",
            filename, package="queeems", mustWork=TRUE)
    }
    return(filePath)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

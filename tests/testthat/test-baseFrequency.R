# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("baseFrequency functions work as intended", {
    sqdatum <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
        "CCTCAGATCACTCTTGTCTCAATAAAAGTA","TGG-AGCGACCCCTTGTCTCAATAAAAATA")
    sqdata <- c("FLDGADKAQECT", "F-EGIDKAQEAT", "FLDGIDYAQEAG")
    aaData <- Biostrings::AAStringSet(sqdata)
    dnaData <- Biostrings::DNAStringSet(sqdatum)
    codonData <- Biostrings::BStringSet(sqdatum)
    names(codonData) <- names(dnaData) <- names(aaData) <- c("s1","s2","s3")
    testaa <- file.path(tempdir(), "testAmino.fasta")
    testdna <- file.path(tempdir(), "testDna.fasta")
    Biostrings::writeXStringSet(dnaData, testdna)
    Biostrings::writeXStringSet(aaData, testaa)
    trial1 <- baseFrequency.noncodon(dnaData)
    trial2 <- baseFrequency.noncodon(aaData)
    trial3 <- baseFrequency.codon(codonData)
    trial4 <- bstringCodons(codonData)
    trial5 <- codonCensusVector( trial4[,1])
    expect_true( max( colSums(trial3)) <= nrow(trial4))
    expect_equal( as.numeric( table( rowSums(trial3))), c(47,2,10,1,1))
    expect_true( all( c(trial2["F",1],trial2["G",4]) == 3))
    expect_true( all( trial5[1,c("CCT","TGG")] == c(2,1)) )
    expect_true( all( trial1[,25] == trial1[,26]))
    expect_identical( trial1[,14], trial1[,29] )
    expect_error( baseCheck("amino") )
    expect_equal( baseCheck("AA"), "aa" )
    expect_identical( baseCheck("DNA"), "dna" )
    expect_identical( baseCheck("CoDoN"), "codon")
    expect_warning( baseFrequency(testdna, "aa") )
    expect_warning( baseFrequency(testaa, "codon") )
    expect_no_warning( baseFrequency(testaa, "aa") )
    expect_no_warning( baseFrequency(testdna, "dna") )
    expect_no_warning( baseFrequency(testdna, "codon") )
    expect_warning( expect_warning( baseFrequency(testaa, "dna") ) )
})


# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

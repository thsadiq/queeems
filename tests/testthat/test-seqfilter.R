# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("seqfilter function works as intended", {
    dsite0 <- sample(seq(20,30), 6)
    dsite1 <- sample(seq(-20,-10), 6)
    expect_error( siteIDcheck(dsite0,5))
    expect_error( siteIDcheck(dsite0,10))
    expect_error( siteIDcheck(dsite1,10))
    expect_equal( siteIDcheck(dsite0,40), dsite0)

    sqdata <- c("FLDGADKAQECT", "F-EGIDKAQEAT", "FLDGIDYAQEAG")
    sdatum <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
                "CCTCAGATCACTCTTGTCTCAATAAAAGTA",
                "TGG-AGCGACCCCTTGTCTCAATAAAAATA")
    aaData <- Biostrings::AAStringSet(sqdata)
    dnaData <- Biostrings::DNAStringSet(sdatum)
    codonData <- Biostrings::BStringSet(sdatum)
    names(codonData) <- names(dnaData) <- names(aaData) <- c("s1","s2","s3")
    testoutput <- file.path(tempdir(), "testfilter.fasta")
    testaa <- file.path(tempdir(), "testAmino.fasta")
    testdna <- file.path(tempdir(), "testDna.fasta")
    Biostrings::writeXStringSet(dnaData, testdna)
    Biostrings::writeXStringSet(aaData, testaa)
    expect_no_error( seqfilter(testaa, c(2,4), "aa"))
    expect_no_error( seqfilter(testdna, sample(30,10), "dna", testoutput))
    expect_equal( width( seqfilter(testdna, sample(10,1), "codon"))[1], 27)
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

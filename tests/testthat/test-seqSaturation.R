# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><>< ================================================================ ><>< #

test_that("all the functions in seqSaturation file works well", {
    testFreqs <- nucPIgen(5, 0.4, 5)
    expect_true( all(nucPIgen(10,0.4,5) >= 0))
    expect_equal( nrow(testFreqs), nrow( nucPIgen(8,0.3,15)))
    expect_true( all( round(colSums( nucPIgen(5,0.4,5)),12) == 1))

    expect_true( siteLLmax(sample(13,4), testFreqs) >= 0)

    sqdatum <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
                 "CCTCAGATCACTCTTGTCTCAATAAAAGTA",
                 "TGG-AGCGACCCCTTGTCTCAATAAAAATA")
    testdna <- file.path(tempdir(), "dnaSQ.fasta")
    dnaData <- Biostrings::DNAStringSet( sqdatum )
    Biostrings::writeXStringSet(dnaData, testdna)
    test1 <- seqsat1(testdna, 5, 0.4, 5, 9, 0.40)
    expect_true( is(test1, "saturateBF"))
    expect_true( length( nonvaries(test1)) >= 0)
    expect_equal( sitecount(test1), length( LogL1(test1)) )
    expect_equal( length( BFs(test1)), length( LogL0(test1)) )

    expect_error(scanapp(1))
    expect_error(scanapp("d"))
    expect_equal(scanapp("a"), scanapp("A"))

    test2a <- seqSaturation(testdna, "A", 5, 0.4, 5)
    expect_equal( length( BFs(test2a)), length( LogL0(test2a)) )
    expect_error(seqSaturation(testdna, "B", 5, 0.4, 5))
    expect_error(seqSaturation(testdna, "C", 5, 0.4, 5))
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

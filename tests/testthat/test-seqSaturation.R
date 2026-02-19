# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><>< ================================================================ ><>< #

test_that("all the functions in seqSaturation file works well", {
    tempbeta <- runif(1)
    fakeBM <- matrix(sample(9,40,TRUE), nrow=4)
    # app1 <- oneBetaLL(tempbeta, fakeBM)
    # expect_true( all(app1 < 0) )
    # expect_equal( length(app1), ncol(fakeBM))

    expect_error( betaspec(20, 0.01, 1))
    expect_error( betaspec(3, 0.01, 10))
    expect_error( betaspec(20, 0.01, 10))
    expect_equal( betaspec(4, 0.5, 1.2), 0)

    expect_equal( length( fubarsplit(20, 0.2, 10)), 20)
    expect_equal( sum( fubarsplit(4, 0.5, 1.2) > 1 ), 2)
    expect_equal( sum( fubarsplit(20, 0.2, 10) <= 1 ), 4)

    # testBSL <- betaSpaceLL(20, 0.7, 10, fakeBM)
    # expect_equal( length(testBSL), ncol(fakeBM))
    # expect_error( betaSpaceLL(20, 0.01, 1, fakeBM))

    sqdatum <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
                 "CCTCAGATCACTCTTGTCTCAATAAAAGTA",
                 "TGG-AGCGACCCCTTGTCTCAATAAAAATA")
    testdna <- file.path(tempdir(), "dnaSQ.fasta")
    dnaData <- Biostrings::DNAStringSet( sqdatum )
    Biostrings::writeXStringSet(dnaData, testdna)
    # test1 <- seqsat1(testdna, 10, 0.2, 10)
    # expect_true( is(test1, "saturateBF"))
    # expect_true( length( nonvaries(test1)) >= 0)
    # expect_equal( sitecount(test1), length( LogL1(test1)) )
    # expect_equal( length( BFs(test1)), length( LogL0(test1)) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

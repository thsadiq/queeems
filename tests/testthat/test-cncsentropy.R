# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><>< ================================================================ ><>< #

test_that("code composed in cncsentropy works", {
    sqdata <- c("CCTCAGATCACTCTTTGGCAACGAAACCTT",
                "CCTCAGATCACTCTTGTCTCAATAAAAGTA",
                "CCC-AGCGACCCCTTGTCTCAATAAAAATA")
    codonData <- Biostrings::BStringSet(sqdata)
    names(codonData) <- c("s1","s2","s3")
    testdna <- file.path(tempdir(), "testDna.fasta")
    Biostrings::writeXStringSet(codonData, testdna)
    test1a <- cncsentropy(testdna, TRUE, 1, "renyi", .1)
    test1b <- cncsentropy(testdna, FALSE, 1, "renyi", .1)
    expect_equal( gentropy(test1a), gentropy(test1b) )
    expect_true( gentropy(test1a) < 0 )

    sqdatm <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
                "CCTCAGATCACTCTTGTCTCAATAAAAGTA",
                "TGG-AGCGACCCCTTGTCTCAATAAAAATA")
    codonData <- Biostrings::BStringSet(sqdatm)
    names(codonData) <- c("s1","s2","s3")
    Biostrings::writeXStringSet(codonData, testdna)
    expect_true( gentropy(cncsentropy(testdna, FALSE, 1, "shannon")) == 0)
    expect_true( gentropy(cncsentropy(testdna, TRUE, 2, "renyi", 0.3)) == 0)

    useq <- queeemsExtdata("RTI.fasta")
    run1 <- cncsentropy(useq, TRUE, 1, "shannon")
    run2 <- cncsentropy(useq, TRUE, 2, "shannon")
    run3 <- cncsentropy(useq, FALSE, 2, "renyi", 0.23)
    run4 <- cncsentropy(useq, FALSE, 3, "renyi", 0.23)
    expect_true(gentropy(run1) <= gentropy(run2))
    expect_true(gentropy(run4) >= gentropy(run3))
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

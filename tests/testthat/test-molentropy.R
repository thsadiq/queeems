# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("calls of the molentropy.R constituent functions work", {
    testFile <- queeemsExtdata("RTI.fasta")
    examFile <- queeemsExtdata("VertCOI_AA.fasta")
    test0 <- molentropy(testFile, "codon", "shannon")
    test1 <- molentropy(examFile, "aa", "renyi", 0.6)
    test2 <- molentropy(testFile, "dna", "renyi", 0.4)
    expect_true( is.na( renyiA(test0)))
    expect_equal( informula(test1), "renyi")
    expect_warning( molentropy(examFile, "codon", "shannon"))
    expect_error( molentropy(testFile, "dna", "renyi", -0.5))
    expect_equal( summary(test1)["noinfo.sites",], length( nonvaries(test1)))
    expect_equal( length(sitentropies(test2)), 3*length(sitentropies(test0)))
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

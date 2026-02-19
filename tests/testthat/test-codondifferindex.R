# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("there are no issues with the codondifferindex function", {
    testFile <- queeemsExtdata("II.fasta")
    test0 <- codondifferindex(testFile, 2, "renyi", 0.4)
    expect_true( renyiA(codondifferindex(testFile, 1, "renyi", 0.4)) == 0.4)
    expect_true( length(sitentropies(test0)) >= length(nonvaries(test0)))
    expect_no_error( codondifferindex(testFile, 1, "renyi", 0.4))
    expect_error( codondifferindex(testFile, 0, "renyi", 0.25))
    expect_no_error( codondifferindex(testFile, 2, "shannon"))
    expect_error( codondifferindex(testFile, 4, "shannon"))
    expect_error( codondifferindex(testFile, 1, "renyi"))
    expect_false( informula( test0) == "shannon")
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

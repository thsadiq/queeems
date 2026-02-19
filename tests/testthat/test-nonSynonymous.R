# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("function in the nonSynonymous file are in order", {
    expect_equal( sum( maxdiff(0)), 61 )
    expect_true( sum( nsynvect("AAA", TRUE)) == sum( nsynvect("AAC", TRUE)))
    expect_equal( nonSynonymous(TRUE, 0), nonSynonymous(FALSE, 0) )
    expect_error( nonSynonymous(TRUE, 4) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

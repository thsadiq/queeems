# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("citing-class.R contents are working fine", {
    test1 <- aboutQueeems()
    expect_no_error( show(test1))
    expect_error( aboutQueeems(9) )
    expect_true( is(test1, "citing"))
    expect_no_error( aboutQueeems(1) )
    expect_error( new("citing", citeText=1))
    expect_no_error( show( new("citing", citeText="Just testing!")) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

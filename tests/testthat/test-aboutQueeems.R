# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("aboutQueeems returns citations as expected", {
    expect_error( aboutQueeems(10) )
    expect_no_error( aboutQueeems() )
    expect_no_error( aboutQueeems(NULL) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

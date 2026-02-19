# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><>< ================================================================ ><>< #

test_that("mnomlogl-related functions respond accordingly", {
    rawpi <- runif(20)
    goodpi <- rawpi / sum(rawpi)
    expect_error( accessinputs(sample(5,8,TRUE), runif(7)) )
    expect_warning( accessinputs(sample(5,8,TRUE), runif(8)) )
    expect_error( accessinputs(sample(5,8,TRUE), runif(8,-1,0)) )
    expect_true( abs(sum(accessinputs(sample(40,20),goodpi)) - 1) < 1e-12)
    expect_true( mnomLogl(sample(40,20),goodpi) < 0)
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

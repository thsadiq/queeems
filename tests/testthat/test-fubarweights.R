# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><>< ================================================================ ><>< #

test_that("all the functions in fubarweights file works well", {
    tempbeta <- runif(1)
    fakeBM <- matrix(sample(9,40,TRUE), nrow=4)

    expect_error( pispecs(20, 0.01, 1))
    expect_error( pispecs(3, 0.01, 10))
    expect_error( pispecs(20, 0.01, 10))
    expect_equal( pispecs(4, 0.5, 1.2), 0)

    expect_equal( fubarweights(30, 0.2, 10, TRUE)[1], 0)
    expect_true( fubarweights(30, 0.2, 10, FALSE)[1] > 0)
    expect_equal( length( fubarweights(20, 0.2, 10, TRUE)), 20)
    expect_equal( length( fubarweights(40, 0.2, 10, FALSE)), 40)
    expect_equal( sum( fubarweights(4, 0.5, 1.2, FALSE) > 1 ), 2)
    expect_equal( sum( fubarweights(20, 0.2, 10, FALSE) <= 1 ), 4)
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

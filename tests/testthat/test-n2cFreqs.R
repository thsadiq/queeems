# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><>< ================================================================ ><>< #

test_that("n2cFreqs.R is in excellent shape", {
    expect_equal( codon2int("CAG"), c(2,1,3))
    expect_equal( codon2int("CGT"), c(2,3,4))

    expect_error( vetFreqs( runif(3)))
    expect_warning( vetFreqs( rep(0.1,5)))
    expect_error( vetFreqs( c(-1, runif(3))))
    expect_error( vetFreqs( runif(4, -5, -1)))
    expect_error( vetFreqs( runif(4, -5, -1)))
    expect_error( vetFreqs( LETTERS[seq(1,4)]))
    expect_error( vetFreqs( c(A=.4,C=.2,G=.2,F=.2)) )
    expect_equal( vetFreqs( c(A=.4,C=.2,G=.2,T=.2)), c(A=.4,C=.2,G=.2,T=.2))

    testrun <- n2cFreqs(c(0.4,0.2,0.2,0.2))
    expect_true(testrun[25] == testrun[10])
    expect_equal( round( sum(testrun), 12), 1)
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

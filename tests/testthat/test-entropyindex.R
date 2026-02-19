# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><>< ================================================================ ><>< #

test_that("code composed in entropyindex function as expected", {
    egcounts <- sample(seq(1,50), 7)
    expect_error( validBases( c(NA,NA,NaN), FALSE))
    expect_error( validBases( c(egcounts,Inf), TRUE))
    expect_no_warning( validBases( c(egcounts,0), TRUE))
    expect_equal( validBases( c(0,egcounts,0), TRUE), egcounts)
    expect_equal( validBases( c(0,egcounts,0), FALSE), c(0,egcounts,0))

    expect_error( siFreq(egcounts,sum(egcounts)-1,TRUE))
    expect_true( sum(siFreq(egcounts,sum(egcounts)+9,TRUE)) < 1)
    expect_true( all( siFreq(egcounts,sum(egcounts),TRUE) <= 1))
    expect_true( all( siFreq(egcounts,sum(egcounts),FALSE) >= 0))
    expect_true( round( sum( siFreq(egcounts,sum(egcounts),TRUE)), 4) == 1)
    expect_true( round( sum( siFreq(egcounts,sum(egcounts),FALSE)), 4) == 1)

    expect_equal( si.shannon(rep(1,4),4), -log(0.25))
    expect_error( si.renyi(egcounts, sum(egcounts), 1.4) )
    expect_error( si.renyi(egcounts, sum(egcounts), -1.4) )

    expect_error( checkFml(NA))
    expect_error( checkFml("logl"))
    expect_error( checkFml("entropy"))
    expect_no_error( checkFml("Renyi"))

    expect_equal( entropyindex(rep(0,4), 10), 0)
    expect_no_error( entropyindex(rep(1,4), 4, ralpha=1))
    expect_equal( entropyindex(rep(0,7), 10, "renyi", 0.2), 0)
    expect_error( entropyindex(rep(1,4), 4, "renyi", ralpha=1))
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

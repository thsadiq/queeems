# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><>< ================================================================ ><>< #

test_that("code compiled in saturateBF-class.R is operational", {
    nFail <- 7
    nFixed <- 4
    nsiteIN <- 12
    h0Total1 <- 10
    h0Total2 <- 1e+04
    h0Total3 <- 1e+05
    altTotal <- 50000
    hsplit <- c( gtools::rdirichlet(1, rexp(nsiteIN,1)))
    bfIN <- c( runif(nFail), runif(nsiteIN-nFail,101,150))
    trialdef1 <- new("saturateBF", bf.pv=bfIN, nuLL=h0Total1*hsplit,
        altLL=altTotal*hsplit, noInfo=sample(nsiteIN,nFixed),
        nsites=nsiteIN, bfcutoff=50, sspi=0.10)
    trialdef2 <- new("saturateBF", bf.pv=bfIN, nuLL=h0Total2*hsplit,
        altLL=altTotal*hsplit, noInfo=sample(nsiteIN,nFixed),
        nsites=nsiteIN, bfcutoff=100, sspi=0.30)
    trialdef3 <- new("saturateBF", bf.pv=bfIN, nuLL=h0Total3*hsplit,
        altLL=altTotal*hsplit, noInfo=sample(nsiteIN,nFixed),
        nsites=nsiteIN, bfcutoff=100, sspi=0.30)
    expect_no_error( show(trialdef1))
    expect_true( all( BFs(trialdef1) > 0))
    expect_equal( nrow( summary(trialdef2)), 5)
    expect_equal( sitecount(trialdef1), nsiteIN)
    expect_equal( sum( LogL0(trialdef2)), h0Total2)
    expect_equal( sum( LogL1(trialdef3)), altTotal)
    expect_equal( length( nonvaries(trialdef1)), nFixed)
    expect_equal( summary(trialdef1)["seq.state",], "Saturated")
    expect_equal( summary(trialdef3)["seq.state",], "Not Saturated")
    expect_equal( summary(trialdef3)["satur.s",], as.character(nsiteIN-nFail))
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

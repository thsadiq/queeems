# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><>< ================================================================ ><>< #

test_that("code compiled in saturatePV-class.R is operational", {
    nPass <- 7
    nFixed <- 4
    nsiteIN <- 12
    h0Total1 <- 10
    h0Total2 <- 1e+05
    altTotal <- 50000
    pvIN <- c( runif(7, 0.05, 1), runif(nsiteIN-nPass,0,0.05))
    hsplit <- c( gtools::rdirichlet(1, rexp(nsiteIN,1)))
    trialdef1 <- new("saturatePV", bf.pv=pvIN, nuLL=h0Total1*hsplit,
        altLL=altTotal*hsplit, noInfo=sample(nsiteIN,nFixed), nsites=nsiteIN)
    trialdef2 <- new("saturatePV", bf.pv=pvIN, nuLL=h0Total2*hsplit,
        altLL=altTotal*hsplit, noInfo=sample(nsiteIN,nFixed), nsites=nsiteIN)
    expect_no_error( show(trialdef1))
    expect_true( all( PVs(trialdef1) >= 0))
    expect_true( all( PVs(trialdef1) <= 1))
    expect_equal( nrow( summary(trialdef2)), 9)
    expect_equal( sitecount(trialdef1), nsiteIN)
    expect_equal( sum( LogL0(trialdef2)), h0Total2)
    expect_equal( sum( LogL1(trialdef2)), altTotal)
    expect_equal( length( nonvaries(trialdef1)), nFixed)
    expect_equal( summary(trialdef2)["seq.state",], "Saturated")
    expect_true( as.numeric( summary(trialdef2)["p.value",]) >= 0.05)
    expect_equal( summary(trialdef1)["satur.s",], as.character(nPass))
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

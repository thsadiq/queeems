# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("none of the code in geneindex-class.R fails", {
    expect_no_error( new("geneindex", eindex=runif(1), noinfo=sample(10,3),
        meth="ReNyI", alphaR=0.2, nbases=10, synonym=TRUE) )
    expect_no_error( new("geneindex", eindex=runif(1), noinfo=sample(10,3),
        meth="ReNyI", alphaR=0.2, nbases=10, synonym=FALSE) )
    expect_no_error( new("geneindex",
        eindex=runif(1), noinfo=seq(1,4)[seq(1,4)>4], meth="SHANNON") )
    testdata <- queeemsExtdata("PI.fasta")
    testrun1 <- cncsentropy(testdata, TRUE, 1, "shannon")
    testrun2 <- cncsentropy(testdata, FALSE, 2, "renyi", 0.23)
    expect_no_error( show(testrun1) )
    expect_false( useSyn(testrun1) == useSyn(testrun2))
    expect_true( gentropy(testrun1) != gentropy(testrun2))
    expect_false( informula(testrun1) == informula(testrun2))
    expect_true( is.na(renyiA(testrun1)) & renyiA(testrun2)==0.23)
    expect_true( sitecount(testrun1) >= length(nonvaries(testrun2)))
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("none of the code in siteindices-class.R fails", {
    simpleEG <- new("siteindices", eindices=runif(10),
        noinfo=c(3,5,7), meth="ReNyI", alphaR=0.2)
    expect_no_error( show(simpleEG))
    expect_no_error( renyiA(simpleEG))
    expect_no_error( summary(simpleEG))
    expect_no_error( informula(simpleEG))
    expect_equal( nonvaries(simpleEG), c(3,5,7))
    expect_equal( length( sitentropies(simpleEG)), 10)
    expect_error( new("siteindices", noname="undefined_entry"))
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

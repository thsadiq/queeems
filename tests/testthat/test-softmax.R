# ><>< ================================================================ ><>< #
# ><><        Saturation Threshold of Codon Evolutionary Models.        ><>< #
# ><>< ================================================================ ><>< #

test_that("softmax transformation produces desired outcomes", {
    testcase <- sample(50,30)
    expect_error( wcheck(rep(0,10)) )
    expect_error( wcheck(rep(NA,10)) )
    expect_equal( wcheck(c(5,1,2,0)), c(5,1,2,0))
    expect_warning( wcheck(c(rep(NA,2), sample(20,18))) )
    expect_equal( round( sum( softmax(testcase, TRUE)), 8), 1.0)
    expect_equal( round( sum( softmax(testcase, FALSE)), 8), 1.0)
    expect_false( all( softmax(testcase, TRUE) == softmax(testcase, FALSE)))
    expect_equal( softmax(sample(10,1), TRUE), softmax(sample(20,1), FALSE))
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

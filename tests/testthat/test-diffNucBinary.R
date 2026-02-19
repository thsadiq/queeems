# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("diffNucBinary-related functions work appropriately", {
    test1 <- diffNucleotide("AAA")
    expect_equal( sum(test1 >= 4), 0 )
    expect_equal( sum(test1 == 0), 1 )
    expect_equal( sum(test1 == 1), 8 )
    expect_equal( sum(test1 == 2), 25 )
    expect_equal( sum(test1 == 3), 27 )
    randID <- sample(61, 1)
    diffCount <- sample(3, 1)
    test2 <- diffNucBinary(diffCount)
    expect_equal( test2[randID,], test2[,randID] )
    expect_equal( sum( diffNucBinary(0)), sum( diag( diffNucBinary(0))) )
    expect_error( diffNucBinary(4) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

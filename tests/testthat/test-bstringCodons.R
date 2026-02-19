# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("bstringCodons-related functions work correctly", {
sqdatum <- c("CCTCAGATCACTCTTTGGCAACGACCCCTT",
    "CCTCAGATCACTCTTGTCTCAATAAAAGTA","TGG-AGCGACCCCTTGTCTCAATAAAAATA")
    codonData <- Biostrings::BStringSet(sqdatum)
    trial <- bstringCodons(codonData)
    expect_equal( dim(trial), c(3,10) )
    expect_warning( codonExtant(LETTERS[1:10]) )
    expect_equal( trial[,3], c("ATC","ATC","CGA") )
    expect_equal( miniGEN(1, c("A","B","C")), "ABC" )
    expect_error( bstringCodons( matrix(LETTERS[1:18],3)) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

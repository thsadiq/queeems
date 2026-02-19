# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("queeemsExtdata functions as expected", {
    expect_true( "RTI_AA.fasta" %in% queeemsExtdata() )
    expect_true( "II.fasta" %in% queeemsExtdata() )
    expect_error( queeemsExtdata("NO_FILE") )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

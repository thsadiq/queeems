# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("codonDissimilarity works as expected", {
    fastaPath <- queeemsExtdata("PI.fasta")
    bstringdatum <- Biostrings::readBStringSet(fastaPath)
    nCodonSites <- unique( width( bstringdatum)) / 3
    expect_equal( nrow( codonDissimilarity(bstringdatum, 1)), 61 )
    expect_equal( ncol( codonDissimilarity(bstringdatum, 2)), nCodonSites )
    expect_error( codonDissimilarity(matrix(NA, 2, 4), 3) )
    expect_error( codonDissimilarity(bstringdatum, 4) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

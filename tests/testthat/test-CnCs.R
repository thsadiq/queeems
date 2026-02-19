# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("all is well with CnCs-related functions", {
    testFile <- queeemsExtdata("II.fasta")
    testdata <- Biostrings::readBStringSet(testFile)
    alignment <- bstringCodons(testdata)
    novaryings <- which( apply(alignment, 2, getinvary))
    rsample <- sample(novaryings, 2)
    case1 <- names( table(alignment[,rsample[1]]))
    wild1 <- case1[length(case1)]
    case2 <- names( table(alignment[,rsample[2]]))
    wild2 <- case2[length(case2)]
    a1 <- nscodoncount(wild1, alignment[,rsample[1]], nonSynonymous(TRUE,0))
    a2 <- nscodoncount(wild2, alignment[,rsample[2]], nonSynonymous(FALSE,1))
    a3 <- nscodoncount(wild1, alignment[,rsample[1]], nonSynonymous(TRUE,2))
    a4 <- nscodoncount(wild2, alignment[,rsample[2]], nonSynonymous(FALSE,3))
    expect_true( all( c(a1,a2,a3,a4) == 0 ))
    expect_true( getinvary(c("___","AAT","AAT")))
    expect_false( getinvary(c("AAC","AAT","AAT")))
    a5 <- nsfreqs( CnCs(testdata, TRUE, 1))[novaryings]
    a6 <- nsfreqs( CnCs(testdata, FALSE, 1))[novaryings]
    expect_true( all(a5 == a6 ) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

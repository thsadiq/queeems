# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("operation of cncsframe-class.R is smooth", {
    tempseq <- queeemsExtdata("VertCOI.fasta")
    moleset <- Biostrings::readBStringSet(tempseq)
    execute <- CnCs(moleset, FALSE, 3)
    execwild <- wildtriplets(execute)
    expect_no_error( show(execute))
    expect_error( CnCs(moleset, TRUE, 4))
    expect_false( nucbalance(execute) < 3 )
    expect_no_error( CnCs(moleset, TRUE, 1))
    expect_false( nORs(execute) == "Synonymous")
    expect_equal( nORs(CnCs(moleset,TRUE,1)), "Synonymous")
    expect_true( all( nsfreqs( CnCs(moleset, FALSE, 0)) == 0) )
    expect_true( length(nonvaries(execute)) <= length(execwild))
    expect_true( all(execwild == wildtriplets( CnCs( moleset, TRUE, 2))))
    expect_false( all( nsfreqs(execute) == nsfreqs(CnCs(moleset,TRUE,2))))

    test2 <- new("cncsframe", nscensus=sample(20,13),
        wildseq=matrix(sample(senseCodon, 13, TRUE),nrow=1),
        nstype="Non-Synonymous", maxdiff=2, invariants=sample(20,2))
    expect_equal( nORs(test2), "Non-Synonymous")
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

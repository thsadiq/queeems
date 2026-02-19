# ><>< ================================================================ ><>< #
# >< Quantify the Extent of Evolutionary Evidence in Molecular Sequences. >< #
# ><>< ================================================================ ><>< #

test_that("none of the code in baseSummary-class.R misbehaves", {
    testfile <- queeemsExtdata("VertCOI.fasta")
    testexe1 <- baseFrequency(testfile,"codon")
    testexe2 <- baseFrequency(testfile, "dna")
    expect_no_error( show(testexe2) )
    expect_warning( baseFrequency(testfile, "aa"))
    expect_equal( nseqs(testexe1), nseqs(testexe2) )
    expect_equal( datatype(testexe1), c(codon="Codon") )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #

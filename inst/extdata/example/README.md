# Folder Information
This folder contains files that were used to generate and process the data sets that are presented in the `queeems` package's vignettes and the manuscript prepared for the Journal of Computational Biology.

- [stemL](stemL): contains the resources for the changing stem length analyses.
 - [seqs03](stemL/seqs03.txt): five replicates of the (100 codon sites)-by-(64 leaves) genetic sequences simulated with each stem length of the balanced evolutionary tree set equal to 0.03.
 - [seqs06](stemL/seqs06.txt): five replicates of the (100 codon sites)-by-(64 leaves) genetic sequences simulated with each stem length of the balanced evolutionary tree set equal to 0.06.
 - [seqs09](stemL/seqs09.txt): five replicates of the (100 codon sites)-by-(64 leaves) genetic sequences simulated with each stem length of the balanced evolutionary tree set equal to 0.09.
 - [seqs12](stemL/seqs12.txt): five replicates of the (100 codon sites)-by-(64 leaves) genetic sequences simulated with each stem length of the balanced evolutionary tree set equal to 0.12.
 - [seqs12](stemL/seqs15.txt): five replicates of the (100 codon sites)-by-(64 leaves) genetic sequences simulated with each stem length of the balanced evolutionary tree set equal to 0.15.
 - [simulate_stem](stemL/simulate_stem.R): R code used to execute the `scoup` simulations that generated the genetic sequences for the changing stem length analyses that are saved in the [stemL](stemL) folder.
 - [summarise_stem](stemL/summarise_stem.R): R code used to perform the `queeems` analyses with respect to changing branch lengths that are reported in the vignettes and the manuscript.

- [varNS](varNS): contains the resources for the changing natural selection effect analyses.
 - [gdata1](varNS/gdata1.txt): five replicates of the (100 codon sites)-by-(64 leaves) genetic sequences simulated with the variance of the non-synonymous selection coefficients set equal to 0.1.
 - [gdata2](varNS/gdata2.txt): five replicates of the (100 codon sites)-by-(64 leaves) genetic sequences simulated with the variance of the non-synonymous selection coefficients set equal to 0.2.
 - [gdata3](varNS/gdata3.txt): five replicates of the (100 codon sites)-by-(64 leaves) genetic sequences simulated with the variance of the non-synonymous selection coefficients set equal to 0.3.
 - [gdata4](varNS/gdata4.txt): five replicates of the (100 codon sites)-by-(64 leaves) genetic sequences simulated with the variance of the non-synonymous selection coefficients set equal to 0.4.
 - [gdata5](varNS/gdata5.txt): five replicates of the (100 codon sites)-by-(64 leaves) genetic sequences simulated with the variance of the non-synonymous selection coefficients set equal to 0.5.
 - [dnds](varNS/dnds.csv): analytical estimates of the magnitude of the active natural selection pressure for all the simulated data replicates and all the variance values considered.
 - [simulate_nvar](varNS/simulate_nvar.R): R code used to execute the `scoup` simulations that generated the genetic sequences saved in the [varNS](varNS) folder.
 - [summarise_nvar](varNS/summarise_nvar.R): R code used to perform the `queeems` analyses that are reported in the vignettes and the manuscript with respect to changing magnitude of positive natural selection effect.


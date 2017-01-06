# introgression_AWT
Almost all the scripts for the introgression AWT paper.

The scripts for the exome capture part (assembly, annotation, etc) is split over:
* https://github.com/singhal/exomeCapture
* https://github.com/CGRL-QB3-UCBerkeley/denovoTargetCapturePopGen

## Exome Capture Scripts
### Pipeline 1: https://github.com/CGRL-QB3-UCBerkeley/denovoTargetCapturePopGen
- `2-ScrubReads`: used to scrub reads
- `3-GenerateAssemblies`: used to assemble reads
- `4-FinalAssembly`: used to assemble across assemblies
### Pipeline 2: https://github.com/singhal/exomeCapture
- `7annotateContigs.pl`: used to match assembles to original targets
- `8initialSNPset.pl`: used to map reads and call SNPs in transcriptome data
- `9clineSNPset.pl`: used to map reads and call SNPs in clinal populations
- `10benchmarking.pl`: used to get data for loci used in benchmarking.
- `11clineAFfiles.pl`: used to get allele frequency data from the cline SNP set 

## Probe Design Scripts
- `divPoly.pl`: calculate divergence & polymorphism for transcriptome data to identify highly-differentiated transcripts
- `dnds.pl`: calculate a (rough) version of dn/ds for transcriptome data to identify transcripts putatively under positive selection
- `exomeCapture.pl`: the script that actually does most of the work; identifies possible exons, filters for GC content etc.
- `finalExons_species.pl`: picks some exons and adds in the other targeted sequence (i.e., mtDNA, UTRs, etc)
- `finalFile_makeShorter.pl`: makes the exon capture array shorter because the original files were overshoots
- `fst.pl`: calculate a (rough) version of Fst for transcriptome data to identify highly-differentiated transcripts
- `initialExonFile.pl`: generated the exon list to be used for downstream selection
- `utrExtract.pl`: extracts UTR sequences to be used in benchtruthing

## Analysis Scripts
- `aliquotDNA.pl`: figured out amounts of DNA to spike into tube to get pooled libraries
- `averageExonMtCoverage.R`: calculated average coverage across each population in each lineage-pair
- `benchtruthing_af.R`: determines how much variance there is in mtDNA SNPs
- `cline_fitting.py`: fits clines to variants; takes in allele frequencies at each variant in each population
- `divergence_contig.py`: estimates divergence & polymorphism for transcriptome data on a per contig (transcript) basis
- `divergence_full.py`: estimates divergence & polymorphism for transcriptome data on a per exon basis
- `dnds.py`: estimates dn/ds for each transcript using PAML
- `fst_contig.py`: estimates Fst for transcriptome data on a per contig (transcript) basis
- `fst_full.py`: estimates Fst for transcriptome data on a per exon basis
- `makeLDfiles.pl`: looks at how cline centers and widths change over genomic space -- allows Anolis to be used as a reference
- `moransI_bootstrap.py`: takes in the files generated by `makeLDfiles.pl` and generates Moran's I estimates as well as doing some bootstrapping
- `pooled_discrepancy.R`: runs the simple R simulations to see how pooling might affect our ability to infer allele frequency accurately

## Figure Scripts
- `paper_figures.ipynb`: includes Python code for most of the figures in the manuscript and some figures that didn't survive
- `correlationCoverage.R`: produces Figure S5 in Supplementary Information
- `coverageAcrossContacts.R`: produces Figure S6 in Supplementary Information
- `cumulativeCoveragePlot.R`: produces Figure S7 in Supplementary Information
- `divpolyCoverage.R`: produces Figure S8 in Supplementary Information
- `gcContentCoverage.R`: produces Figure S9 in Supplementary Information
- `specificity.R`: produces Figure S4 in Supplementary Information

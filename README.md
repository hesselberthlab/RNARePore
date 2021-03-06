# RNARePore
Ligation by plant and fungal RNA ligases yields an internal 2′-phosphate group on each RNA ligation product. Here, we define several unique signals produced by 2′-phosphorylated RNAs during nanopore sequencing:
* A 2′-phosphate at the splice junction of HAC1 mRNA inhibits 5′→3′ degradation, enabling detection of decay intermediates by 5′ end mapping in yeast RNA repair mutants.
* During direct RNA sequencing, intact 2′-phosphorylated RNAs produce diagnostic changes in nanopore current properties and base calling features.

This repository describes our general bioinformatic strategy for detecting 2′-phosphate signals, as well as the R markdown files used to generate the figures in the preprint.

Alignment references used in this work were:
* An [S. cerevisiae transcriptome reference generated by the Jacobson Lab](https://github.com/Jacobson-Lab/yeast_transcriptome_v5), which contains both pre- and post-spliced mRNA references
* A custom _S. cerevisiae_ tRNA reference containing a consensus sequence for each cytoplasmic tRNA species in budding yeast, built from tRNAscan-SE gene predictions available at [gtRNAdb](http://gtrnadb.ucsc.edu). This reference is specific to the sequencing library preparation described in in [PMID: 34618430](https://pubmed.ncbi.nlm.nih.gov/34618430/), with tRNA splint adapter sequences appended and prepended to enable alignment to both mature tRNA sequence and ligated adapters. All tRNAs in this reference have had CCA added to their 3′ ends, and (where applicable) introns removed in silico.
* 18S and 25S rRNA fastas for reanalysis of dwell time at select sites in [PMID: 35252946](https://pubmed.ncbi.nlm.nih.gov/35252946/) were obtained from the authors' [RNA_nanoSHAPE](https://github.com/physnano/rRNA_nanoSHAPE) repository
* The RNA oligonucleotide sequences in [PMID: 34893601](https://pubmed.ncbi.nlm.nih.gov/34893601/) share a common sequence, but are decorated with different collections of RNA modifications. Their sequence was extracted from the manuscript and can be found here as **oligos.fa**.
* A reference for the synthetic RNA oligonucleotide generated by the splint ligation strategy in Figure S4 is contained in **splint.fa**

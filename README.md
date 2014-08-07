ExomeSeqPipeline
================
Collection of utilities used for secondary and tertiary analysis of Exome Seq data.

## Directory List
* bin
  * jacquard.py
  * jacquard_utils.py
  * pivot_variants.py : consolidates a set of VCF files into a single variant-level Excel file
  * rollup_genes.py : summarizes a variant-level Excel file to create a gene-level Excel file
  * style.py : styles Excel file
  * tag.py : adds Jacquard-specific FORMAT tags to each VarScan or MuTect VCF record
  * normalize_varscan.py : concatenates VarScan snp and indel VCF files into a merged VCF file, adding high-confidence labels
  * normalize_strelka.py : concatenates Strelka snvs and indels VCF files into a merged VCF file
  * normalize_utils.py
  * filter_somatic.py : filters VCF files to contain only variant calls for coordinates called as somatic in any sample
  * consensus.py : calculates consensus calls from VCf file
  * merge.py : merges several Jacquard-tagged VCF files into one VCf file
* test
  * automated tests
  
================

Jessica Bene, Chris Gates - UM BRCF Bioinformatics Core


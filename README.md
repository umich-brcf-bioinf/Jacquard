ExomeSeqPipeline
================
Collection of utilities used for secondary and tertiary analysis of Exome Seq data.

## Directory List
* bin
  * jacquard.py
  * pivot_variants.py : consolidates a set of VCF files into a single variant-level Excel file
  * rollup_genes.py : summarizes a variant-level Excel file to create a gene-level Excel file
  * style.py : stylizes Excel file
  * tag.py : adds Jacquard-specific FORMAT tags to each VarScan or MuTect VCF record
  * normalize_varscan: concatenates VarScan snp and indel VCF files into a merged VCF file
* test
  * automated tests
  
================

Jessica Bene, Chris Gates - UM BRCF Bioinformatics Core


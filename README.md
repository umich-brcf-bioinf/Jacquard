
ExomeSeqPipeline
================
Collection of utilities used for secondary and tertiary analysis of Exome Seq data.

[![Build Status](https://travis-ci.org/umich-brcf-bioinf/ExomeSeqPipeline.svg?branch=develop)](https://travis-ci.org/umich-brcf-bioinf/ExomeSeqPipeline)
[![Coverage Status](https://coveralls.io/repos/umich-brcf-bioinf/ExomeSeqPipeline/badge.png?branch=develop)](https://coveralls.io/r/umich-brcf-bioinf/ExomeSeqPipeline?branch=develop)

## Directory List
* bin
  * jacquard.py
  * jacquard_utils.py
  * normalize.py : concatenates VarScan snp/indel VCF files or Strelka snvs/indels VCF files into a merged VCF file, adding high-confidence labels to VarScan files
  * tag.py : adds Jacquard-specific FORMAT tags to each VarScan or MuTect VCF record
  * filter_somatic.py : filters VCF files to contain only variant calls for coordinates called as somatic in any sample
  * merge.py : merges several Jacquard-tagged VCF files into one VCf file
  * consensus.py : calculates consensus calls from VCf file
  * rollup_genes.py : summarizes a variant-level Excel file to create a gene-level Excel file
  * style.py : styles Excel file
  * pivot_variants.py : consolidates a set of VCF files into a single variant-level Excel file
* test
  * automated tests
  
================

Jessica Bene, Chris Gates - UM BRCF Bioinformatics Core


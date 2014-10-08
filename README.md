
Jacquard
================
Suite of command-line tools to expedite analysis of exome variant data from multiple patients and multiple variant callers.

[![Build Status](https://travis-ci.org/umich-brcf-bioinf/Jacquard.svg?branch=develop)](https://travis-ci.org/umich-brcf-bioinf/Jacquard) 
[![Coverage Status](https://img.shields.io/coveralls/umich-brcf-bioinf/Jacquard.svg)](https://coveralls.io/r/umich-brcf-bioinf/Jacquard)
## Directory List
* jacquard-runner.py : Convenience wrapper for running Jacquard directly from source tree.
* jacquard : Python libraries 
* test : Automated unit tests

## Usage

jacquard-runner.py <subcommand> [options] [arguments]

For help on a specific subcommand:
jacquard-runner.py <subcommand> --help 

## Subcommands
    normalize           Accepts a single directory containing:
                            MuTect VCFs
                            Strelka snvs/indels VCFs
                            VarScan snp/indel VCFs and filtered native output
                        Creates a new directory of merged, sorted VCFs (adding 
                        high confidence tags for merged VarScan results).
    tag                 Accepts a directory of VCF results and creates a new
                        directory of VCFs, adding Jacquard-specific FORMAT
                        tags for each VCF record.
    filter_hc_somatic   Accepts a directory of Jacquard-tagged VCF results
                        from one or more callers and creates a new directory
                        of VCFs, where rows have been filtered to contain only
                        positions that were called as high-confidence somatic 
                        in any VCF.
    merge               Accepts a directory of VCFs and returns a single
                        merged VCF file.
    consensus           Accepts a Jacquard-merged VCF file and creates a new VCF
                        file, adding consensus fields.
    expand              Transforms VCF file into tab-separated text file 
                        expanding INFO fields and FORMAT tags into distinct
                        columns.

================

Jessica Bene, Ashwini Bhasi, Chris Gates, Kevin Meng, Peter Ulintz - UM BRCF Bioinformatics Core


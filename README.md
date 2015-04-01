
#Jacquard
________
Suite of command-line tools to expedite analysis of exome variant data from multiple patients and multiple variant callers.

[![Build Status](https://travis-ci.org/umich-brcf-bioinf/Jacquard.svg?branch=develop)](https://travis-ci.org/umich-brcf-bioinf/Jacquard) 
[![Coverage Status](https://coveralls.io/repos/umich-brcf-bioinf/Jacquard/badge.png?branch=develop)](https://coveralls.io/r/umich-brcf-bioinf/Jacquard?branch=develop)

### Directory List
* jacquard-runner.py : Convenience wrapper for running Jacquard directly from source tree.
* jacquard : Python libraries
* spikes : Unsupported prototypes and other experiements
* test : Automated unit tests

### Usage

jacquard <subcommand> [options] [arguments]

For help on a specific subcommand:
jacquard-runner.py <subcommand> --help 

#### Subcommands
    translate           Accepts a directory of VCF results (including VarScan high
                        confidence files). Creates a new directory of VCFs,
                        adding Jacquard-specific FORMAT tags for each VCF
                        record.
    filter_hc_somatic   Accepts a directory of Jacquard-tagged VCF results
                        from one or more callers and creates a new directory
                        of VCFs, where rows have been filtered to contain only
                        positions that were called as high-confidence somatic 
                        in any VCF.
    merge               Accepts a directory of VCFs and returns a single
                        merged VCF file.
    summarize           Accepts a Jacquard-merged VCF file and creates a new VCF
                        file, adding summary fields/tags.
    expand              Transforms VCF file into tab-separated text file 
                        expanding INFO fields and FORMAT tags into distinct
                        columns.

================

Email bfx-jacquard@umich.edu for support and questions.

UM BRCF Bioinformatics Core



Jacquard
========
Suite of command-line tools to expedite analysis of exome variant data from multiple patients and multiple variant callers.

.. image:: https://travis-ci.org/umich-brcf-bioinf/Jacquard.svg?branch=develop
    :target: https://travis-ci.org/umich-brcf-bioinf/Jacquard
    :alt: Build Status

.. image:: https://coveralls.io/repos/umich-brcf-bioinf/Jacquard/badge.png?branch=develop
    :target: https://coveralls.io/r/umich-brcf-bioinf/Jacquard?branch=develop
    :alt: Coverage Status
    
.. image:: https://pypip.in/license/jacquard/badge.png
    :target: https://pypi.python.org/pypi/jacquard/
    :alt: License

.. image:: https://pypip.in/d/jacquard/badge.png
    :alt: Downloads Counter

Files
-----
* jacquard-runner.py : Convenience wrapper for running Jacquard directly from source tree.
* jacquard : Python libraries
* spikes : Unsupported prototypes and other experiements
* test : Automated unit tests

Usage
-----
``$jacquard <subcommand> [options] [arguments]``

*Subcommands*

:translate:
   Accepts a directory of VCF results (including VarScan high confidence 
   files). Creates a new directory of VCFs, adding Jacquard-specific FORMAT 
   tags for each VCF record.
:merge:
   Accepts a directory of VCFs and returns a single merged VCF file.
   Optionally filters to a subset of variants/loci.
:summarize:
   Accepts a Jacquard-merged VCF file and creates a new VCF file, adding 
   summary fields/tags.
:expand:
   Transforms VCF file into tab-separated text file expanding INFO fields and 
   FORMAT tags into distinct columns.

For help on a specific subcommand:

``jacquard <subcommand> --help``


---

Email bfx-jacquard@umich.edu for support and questions.

UM BRCF Bioinformatics Core


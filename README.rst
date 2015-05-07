========
Jacquard
========
Suite of command-line tools to expedite analysis of exome variant data from multiple patients and multiple variant callers.

.. image:: https://travis-ci.org/umich-brcf-bioinf/Jacquard.svg?branch=develop
    :target: https://travis-ci.org/umich-brcf-bioinf/Jacquard
    :alt: Build Status

.. image:: https://coveralls.io/repos/umich-brcf-bioinf/Jacquard/badge.png?branch=develop
    :target: https://coveralls.io/r/umich-brcf-bioinf/Jacquard?branch=develop
    :alt: Coverage Status

.. image:: https://img.shields.io/pypi/l/Jacquard.svg
    :target: https://pypi.python.org/pypi/jacquard/
    :alt: License

.. image:: http://img.shields.io/pypi/v/colour.svg?style=flat
   :target: https://pypi.python.org/pypi/jacquard/
   :alt: Latest PyPI version

.. image:: https://img.shields.io/pypi/dm/Jacquard.svg
   :target: https://pypi.python.org/pypi/jacquard/
    :alt: Downloads Counter

Files
=====
- jacquard-runner.py : Convenience wrapper for running Jacquard directly from source tree.
- jacquard : Python libraries
- spikes : Unsupported prototypes and other experiments
- test : Automated unit tests

Usage
=====
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
   FORMAT tags into discrete columns.

For help on a specific subcommand:

``jacquard <subcommand> --help``


====


Email bfx-jacquard@umich.edu for support and questions.

UM BRCF Bioinformatics Core


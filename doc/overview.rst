Overview
========

Jacquard, a suite of Python command line tools, provides a practical approach to
integrating complex somatic variant data sets.  It is an open source tool for
expediting variant analysis across multiple patient samples and multiple
variant callers.


The tools is designed to be used by bioinformatic analysts (whoever ran the
variant callers); the output is intended to be useful to analysts and
biological researchers.


Most variant callers have embraced the Variant Call Format (VCF) standard [REFERENCE] which
clearly and succinctly describes variants from a single tumor-normal pair.
However, while many callers follow the standard, they often adopt different
ways to partition results (e.g. somatic file vs. germline file, or
SNP vs. indel); likewise each caller creates its own dialect of VCF fields and
tags. Jacquard transforms the dialects of different variant callers into a
controlled vocabulary of tags with a consistent representation of values.
Furthermore it intelligently merges VCFs from different patients and callers to
create a single unified VCFs across your dataset.

The consistent tag names and represntations expedite downstream analysis and 
the ingerated VCF highlights prevelance of specific variants, and overall
mutation loads across samples.

<schematic picture of excel with highlighted columns/rows (from the poster)>

At this time, the Jacquard-supported variant callers are MuTect, VarScan, and
Strelka. A subset of the Jacquard commands support VCFs from other variant
callers.


Email bfx-jacquard@umich.edu for support and questions.

UM BRCF Bioinformatics Core


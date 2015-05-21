Overview
========

Jacquard, a suite of Python command line tools, provides a practical approach to
integrating complex somatic variant data sets. It is an open source tool for
expediting variant analysis across multiple patient samples and multiple
variant callers.


The tools are designed to be used by bioinformatic analysts (whoever ran the
variant callers); the output is intended to be useful to analysts and
biological researchers.


Most variant callers have embraced the Variant Call Format (VCF) standard
[Reference]_, which clearly and succinctly describes variants from a single
tumor-normal pair. However, while many callers follow the standard, they often
adopt different ways to partition results (e.g. somatic file vs. germline file,
or SNP vs. indel); likewise, each caller creates its own dialect of VCF fields
and tags. Jacquard transforms the dialects of different variant callers into a
controlled vocabulary of tags with a consistent representation of values.
Furthermore, it intelligently merges VCFs from different patients and callers
to create a single, unified VCF across your dataset.

The consistent tag names and represntations expedite downstream analysis, and
the ingerated VCF highlights both the prevelance of specific variants and the
overall mutation loads across samples.

.. figure:: images/overview_Diagram.jpg

   **Overview of Jacquard Workflow :** *Jacquard transforms different caller 
   dialects into a uniform VCF format.*

At this time, the Jacquard-supported variant callers are MuTect, VarScan, and
Strelka. A subset of the Jacquard commands support VCFs from other variant
callers.


**Most Recent Version:** Jacquard-0.41

**Application Size:** jacquard-0.41.tar.gz (md5) : 37KB

**Technology Stack:** Python 2.7, 3.x; natsort 3.5.2; numpy 1.7.1;
testfixtures 3.0.2

|
|
Contact Us
----------

Email bfx-jacquard@umich.edu for support and questions.

**UM BRCF Bioinformatics Core**


.. [Reference] Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, et
   al. The variant call format and VCFtools.Bioinformatics 2011; 27: 2156–8.

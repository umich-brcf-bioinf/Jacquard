Command Details
===============

Jacquard is a suite of tools that can be either run in succession or
individually; the typical workflow is to run:

.. toctree::
   :titlesonly:

   translate
   merge
   summarize
   expand


Translate and summarize are useful only for supported callers; merge and expand
work for any VCFs. Each of these commands is described in detail in the
following pages.


**General usage**
::

   jacquard <SUBCOMMAND> [ARGUMENTS] [OPTIONS]

For help on a specific command:
::

   jacquard <SUBCOMMAND> --help


* Jacquard first writes output files to a temporary directory and only copies
  the files upon successful completion of each subcommand.
* Error, warning, and info messages are written to console and log file. Debug
  messages are only written to the log file (unless --verbose specified).


**Input File Conventions**

* Jacquard assumes that the first element of the filename (up to the first dot)
  is a patient identifier. For example:

 * patientA-113.mutect.vcf
 * patientA-113.strelka.snv.vcf
 * patientA-113.strelka.indel.vcf

  This set of three files all have the same patient identifier (patientA-113).
  The tumor-normal sample pairs will be combined into a single pair of 
  tumor-normals columns in the merged VCF. See :ref:`merge <merge-command>` for
  more details.

* To translate a specific VCF dialect, Jacquard determines the source variant
  caller based on the VCF metaheaders. For this reason it is essential that you
  preserve all metaheaders in the source VCF.

* For a specific source VCF, Jacquard automatically determines the tumor and
  normal samples based on the column header and the metaheaders.



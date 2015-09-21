Command Details
===============

Jacquard is a suite of tools that can be either run in succession or
individually: translate, merge, summarize, and expand. Each of these tools is
discussed in detail in subsequent pages:

.. toctree::
   :titlesonly:

   translate
   merge
   summarize
   expand

The typical workflow for Jacquard is to run VCF files through translate, merge,
summarize, and expand. However, any of the below workflows may be implemented
to obtain meaningful results.

.. figure:: images/general_workflows1.jpg
   
   **Multiple Jacquard Workflows :** *There are 8 different possible workflows
   in Jacquard. The first workflow represented in the above diagram is the 
   signature workflow of Jacquard.*


The Jacquard-produced output VCFs are fully compliant VCF files that can be
easily loaded into an external program, such as an annotation tool.


Jacquard first writes output files to a temporary directory and only copies the
files upon successful completion of each command.


Error, warning, and info messages are written to console and log file. Debug
messages are only written to the log file unless logger is initialized as
verbose (in which case debug is also echoed to console). 


General usage
^^^^^^^^^^^^^
``usage: jacquard <SUBCOMMAND> [ARGUMENTS] [OPTIONS]``


For help on a specific command:


``jacquard <SUBCOMMAND> --help``


* Jacquard first writes output files to a temporary directory and only copies
  the files upon successful completion of each subcommand.
* Error, warning, and info messages are written to console and log file. Debug
  messages are only written to the log file (unless --verbose specified).


Input File Conventions
^^^^^^^^^^^^^^^^^^^^^^
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


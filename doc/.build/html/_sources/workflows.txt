Workflows
=========
Jacquard is a suite of tools that can be either run in succession or
individually: translate, merge, summarize, and expand. Each of these tools is
discussed in detail in subsequent pages.

The typical workflow for Jacquard is to run VCF files through translate, merge,
summarize, and expand. However, any of the below workflows may be implemented
to obtain meaningful results.

.. figure:: images/general_workflows1.jpg
   
   **Multiple Jacquard Workflows :** *There are 8 different possible workflows
   in Jacquard. The first workflow represented in the above diagram is the 
   signature workflow of Jacquard.*


Supported Variant Callers
=========================
Merge and Expand are able to process VCF files from any variant caller.
Translate and Summarize, however, must be run with VCF files from one or more
of the supported variant callers. Currently, Jacquard supports MuTect, VarScan,
and Strelka.

.. figure:: images/general_workflows2.jpg

   **Supported Variant Callers for Jacquard :** *Each subcommand in Jacquard 
   currently supports MuTect, VarScan and Strelka. Merge and Expand commands 
   can also accept VCF files from other variant callers.*
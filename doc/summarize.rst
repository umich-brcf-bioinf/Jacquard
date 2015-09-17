.. _summarize-command:

Summarize
=========
The summarize command adds new INFO fields and FORMAT tags that combine variant
data from the merged VCF. It will only work with VCF files that have been
translated.

.. figure:: images/summarize.jpg

   **Summarizing Format Tags :** *The Jacquard-translated format tags from
   each caller are aggregated and processed together to create consensus format
   tags.* 

Usage
-----
``usage: jacquard summarize <input_file> <output_file>``


*positional arguments:*

+--------+---------------------------------------------------------------------+
| input  | | Jacquard-merged VCF file (or any VCF with Jacquard tags; e.g.     |
|        |   JQ_SOM_MT)                                                        |
+--------+---------------------------------------------------------------------+
| output | | A single VCF file                                                 |
+--------+---------------------------------------------------------------------+

Description
-----------
The summarize command uses the Jacquard-specific tags to aggregate caller
information from the file, providing a summary-level view. The inclusion of
summary fields, such as averages, helps you to easily determine which are the
true variants.

The summarized format tags contain the prefix 'JQ_SUMMARY'.

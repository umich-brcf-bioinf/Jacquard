.. _summarize-command:

Summarize
=========
The summarize command adds new INFO fields and FORMAT tags that combine variant
data from the merged VCF. It will only work with VCF files that have been
translated.

.. figure:: images/summarize.jpg

   **Summarizing Format Tags :** *The Jacquard-translated format tags from
   each caller are aggregated to create summary format tags.* 

Usage
-----
::

   jacquard summarize <input_file> <output_file>


*positional arguments:*

+-------------+----------------------------------------------------------------+
| input_file  | | Jacquard-merged VCF file (or any VCF with Jacquard tags; e.g.|
|             |   JQ_SOM_MT)                                                   |
+-------------+----------------------------------------------------------------+
| output_file | | A single VCF file                                            |
+-------------+----------------------------------------------------------------+

Description
-----------
The summarize command uses the Jacquard-specific tags to aggregate caller
information from the file, providing a summary-level view. The inclusion of
summary fields, such as averages, helps you to easily determine which are the
true variants.

The summarized format tags contain the prefix 'JQ_SUMMARY'.

Example summary FORMAT tags
^^^^^^^^^^^^^^^^^^^^^^^^^^^

+-----------------------+------------------------------------------------------+
| Tag name              | Description                                          |
+=======================+======================================================+
| JQ_SUMMARY_AF_AVERAGE | | Average allele frequency across recognized variant |
|                       | | callers that reported frequency for this position  |
|                       | | [average(JQ_*_AF)].                                |
+-----------------------+------------------------------------------------------+
| JQ_SUMMARY_AF_RANGE   | | Max(allele frequency) - min (allele frequency)     |
|                       | | across recognized callers.                         |
+-----------------------+------------------------------------------------------+
| JQ_SUMMARY_HC_GT      | | High confidence consensus genotype (inferred from  |
|                       | | JQ_*_GT and JQ_*_CALLER_PASSED). Majority rules;   |
|                       | | ties go to the least unusual variant (0/1>0/2>1/1).|
|                       | | Variants which failed their filter are ignored.    |
+-----------------------+------------------------------------------------------+
| JQ_SUMMARY_SOM_COUNT  | | Count of recognized variant callers that reported  |
|                       | | confident somatic call for this sample-position.   |
+-----------------------+------------------------------------------------------+

See summary VCF metaheaders for full list of summary tags and descriptions.


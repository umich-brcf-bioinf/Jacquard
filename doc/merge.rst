.. _merge-command:

Merge
=====
The merge command integrates a directory of VCFs into a single VCF. It is
caller-agnostic and can be used on any set of VCF files.

.. figure:: images/merge_join_step.jpg

   **The Merging Process :** *Sample-specific information is grouped together
   for each patient.*

Usage
-----
::

   jacquard merge <input_dir> <output_file>
                  [--include_format_tags=JQ_.*]
                  [--include_cells=valid]
                  [--include_rows=at_least_one_somatic]


*positional arguments:*

+-------------+----------------------------------------------------------------+
| input_dir   | | Directory containing input VCF files ti be merged            |
+-------------+----------------------------------------------------------------+
| output_file | | An integrated VCF file                                       |
+-------------+----------------------------------------------------------------+


*optional arguments:*

+------------------------+-----------------------------------------------------+
| --include_format_tags= | | Comma-separated user-defined list of regular      |
|                        | | expressions for format tags to be included in     |
|                        | | output; (defaults to **'JQ_.*'**)                 |
+------------------------+-----------------------------------------------------+
| --include_cells=       | | all:  Include all variants                        |
|                        | | **valid**:  Only include valid variants           |
|                        | | passed:  Only include variants which passed their |
|                        | |          respective filter                        |
|                        | | somatic:  Only include somatic variants           |
+------------------------+-----------------------------------------------------+
| --include_rows=        | | all:  Include all variants at loci                |
|                        | | at_least_one_passed:  Include all variants at loci|
|                        | |    where at least one variant passed              |
|                        | | all_passed:  Include all variants at loci where   |
|                        | |    all variants passed                            |
|                        | | **at_least_one_somatic**:  Include all variants at|
|                        | |    loci where at least one variant was high-      |
|                        | |    confidence somatic                             |
|                        | | all_somatic:  Include all variants at loci where  |
|                        | |    all variants were high-confidence somatic      |
+------------------------+-----------------------------------------------------+
| --include_all          | | Equivalent to:                                    |
|                        | |    --include_format_tags='.*'                     |
|                        | |    --include_cells=all                            |
|                        | |    --include_rows=all                             |
|                        | | Useful when merging untranslated VCFs which have  |
|                        | | already been filtered to passing variants.        |
+------------------------+-----------------------------------------------------+



Description
-----------
Conceptually, merge has four basic steps, each described in detail below:
 #. Integrate matching loci from different VCFs into common rows
 #. Combine matching samples from different VCFs into common columns
 #. Filter tag values and rows
 #. Assemble the subset of FORMAT tags to be included in the final VCF

Integrate matching loci
^^^^^^^^^^^^^^^^^^^^^^^
Jacquard first develops the superset of all loci (CHROM, POS, REF, and ALT) 
across the set of all input VCFs. For each locus, the input VCF FORMAT tags and
values are merged into a single row. Input variant record-level fields (such as
FILTER, INFO, etc.) are ignored.

MERGE_LOCI_IMAGE_HERE


Combine matching samples
^^^^^^^^^^^^^^^^^^^^^^^^
In the input directory, an individual sample could be called by more than one
variant caller. When merging, Jacquard combines results from the same sample
into a single column. Merged sample names are constructed by concatenating the
filename prefix and the VCF column header.

+--------------------+-----------------------------------+---------------------+
| Filename           | VCF Column header                 | Merged sample names |
+--------------------+-----------------------------------+---------------------+
| case_A.strelka.vcf | #CHROM ... FORMAT SAMPLE1 SAMPLE2 | | case_A:SAMPLE1    |
|                    |                                   | | case_A:SAMPLE2    |
+--------------------+-----------------------------------+---------------------+
| case_A.mutect.vcf  | #CHROM ... FORMAT SAMPLE1 SAMPLE2 | | case_A:SAMPLE1    |
|                    |                                   | | case_A:SAMPLE2    |
+--------------------+-----------------------------------+---------------------+
| case_B.strelka.vcf | #CHROM ... FORMAT SAMPLE3 SAMPLE4 | | case_B:SAMPLE3    |
|                    |                                   | | case_A:SAMPLE4    |
+--------------------+-----------------------------------+---------------------+
| case_B.mutect.vcf  | #CHROM ... FORMAT SAMPLE3 SAMPLE4 | | case_B:SAMPLE3    |
|                    |                                   | | case_A:SAMPLE4    |
+--------------------+-----------------------------------+---------------------+

Given the input VCFs above, the resulting merged VCF will have four sample
columns:
case_A|SAMPLE1,  case_A|SAMPLE2,  case_B|SAMPLE1,  case_B|SAMPLE2.


Filter tag values and rows
^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, merge contains only Jacquard-translated format tags (JQ\_\.*) and
includes all variants with valid syntax at loci where at least one variant was
somatic. The resulting filtered files contain fewer rows, yet higher quality
data than the input files.

Though most variant callers have their own distinct set of format tags, some
tag names are common across multiple callers. If there are any format tag name
collisions, merge will add a prefix (e.g. JQ1_<original_tag>) in order to
disambiguate the format tags.


.. figure:: images/merge_filter_step.jpg

   **The Filtering Process :** *Rows and specific cells in the VCF files are 
   filtered based on the command-line options.*

After filtering, the merge command combines all of the input VCFs into a single,
merged VCF that includes all necessary information for continuing your analysis.

The resulting VCF files contain the distinct set of all coordinates (CHROM, POS,
REF, and ALT) and samples from the input files, provided they pass the filters.
Each coordinate from the input VCF files is added to the output file, which
increases the file length. Additionally, sample columns are merged for each
patient, adding sample specific information and leading to increased column and
file width.

.. note:: Importantly, rather than giving caller-wise sample columns in the
          output VCF file, merge emits patient-wise sample columns. For each
          patient, the merge command joins the set of corresponding sample
          columns into a single column. The grouping of sample-specific
          information for each patient helps to easily analyze the data.


Assemble the subset of FORMAT tags
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO


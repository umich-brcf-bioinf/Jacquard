Translate
---------
The translate command accepts a directory of VCF files and creates a new
directory of "translated" VCF files, which include several Jacquard-specific
format tags and their corresponding metaheaders. Were a variant in the source
VCF to be malformed, the translated VCF file would label it as such in the
FILTER column.


You can gather all input VCFs into a single directory and run translate once, or
partition VCFs into separate directories (for example, by variant caller) and
run translate once for each input directory. When partitioning into separate
input directories, all file names must be unique.


The translated format tags contain a caller specific prefix; example: 'JQ_SK'
for Strelka, 'JQ_VS' for VarScan and 'JQ_MT' for MuTect.

.. figure:: images/translate_pic.jpg

   **Addition of the Jacquard-Specific Format Tags :** *The translated VCF files 
   contain the original format tags from the input files as well as the 
   Jacquard-specific format tags.*

|

Usage
^^^^^
``jacquard translate <input_dir> <output_dir> [OPTIONS]``


*positional arguments:*

=====================================  ========================================
input                                  Directory containing VCF files (and 
                                       VarScan high confidence files). Other
                                       file types ignored
output                                 Directory containing VCF files. Will
                                       create if doesn't exist and will
                                       overwrite files in output directory as
                                       necessary
=====================================  ========================================


*optional arguments:*

=====================================  ========================================
--allow_inconsistent_sample_sets
                                       Set this flag if not every patient is
                                       represented by the same set of
                                       caller-VCFs.
--varscan_hc_filter_file_regex REGEX   Regex pattern that identifies optional
                                       VarScan high-confidence filter files.
                                       The VCF, high-confidence file pairs
                                       should share the same prefix. For
                                       example, given patientA.snp.vcf,
                                       patientA.indel.vcf,
                                       patientA.snp.fpfilter.pass, and
                                       patientA.indel.fpfilter.pass, you could
                                       enable this option as
                                       varscan_hc_filter_file_regex=
                                       '.fpfilter.pass$'
=====================================  ========================================

   VarScan tips:
   To translate VarScan calls, Jacquard requires the VarScan VCF files (snp
   and/or indel). Jacquard can additionally accept VarScan somatic high-
   confidence files; the translate command supresses VCF variant records that
   are absent in the high-confidence files. Somatic high-confidence files
   should be placed alongside corresponding VarScan VCFs and must have the same
   prefix as their corresponding VCF file. The suffix should be specified using
   the command line argument.

   Note: Optional high-confidence files are not VCF files.

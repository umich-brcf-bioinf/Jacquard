General usage
=============
<General data flow diagram>

usage: jacquard <SUBCOMMAND> [ARGUMENTS] [OPTIONS]

For help on a specific command:
 'jacquard <SUBCOMMAND> -h' 

subcommands: {translate, merge, summarize, expand} 

translate   Accepts a directory of VCF files (and VarScan high confidence
            files). Creates a new directory of VCF files, adding
            Jacquard-specific FORMAT tags to each VCF record.

merge       Accepts a directory of VCF files and returns a single merged VCF
            file.

summarize   Accepts a Jacquard-tagged VCF file and creates a new file, adding
            summary fields. 

expand      Converts a VCF file into an tab-delimited file in a tabular format
            which would be better suited than a VCF file for analysis and
            visualization in R, Pandas, Excel, or another third party
            application.

See https://github.com/umich-brcf-bioinf/Jacquard for more information.


Input file conventions
^^^^^^^^^^^^^^^^^^^^^^
Jacquard assumes that the first element of the filename is a patient identifier.

   patientA-113.mutect.vcf
   patientA-113.strelka.snv.vcf
   patientA-113.strelka.indel.vcf

This set of three files all have the same patient identifier (patientA-113) to
represent the same tumor-normal pair; these files will be combined into a
single pair of tumor-normals in the merged VCF. See merge for more details.


To translate a specific VCF dialect, Jacquard determines the source variant
caller based on the VCF metaheaders. For this reason it is essential that you
preserve all metaheaders in the source VCF.


For a specific source VCF, Jacquard automatically determines the tumor and
normal samples based on the column header and the metaheaders.


Notes
^^^^^
The Jacquard-produced output VCFs are fully compliant VCF files that can be
easily loaded into an external program, such as an annotation tool.


Jacquard first writes output files to a temporary directory and only copies the
files upon successful completion of each subcommand.


Error, warning, and info messages are written to console and log file. Debug
messages are only written to the log file unless logger is initialized as
verbose (in which case debug is also echoed to console). 


Sub-commands
============

Common Options
^^^^^^^^^^^^^^
  -h, --help            Show this help message and exit
  --force               Overwrite contents of output directory
  --log_file LOG_FILE   Path to log file destination. Defaults to current 
                        working directory if not specified. 
  -v, --verbose         Prints detailed logs to the console


Translate
^^^^^^^^^
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

   <translate diagram>

Usage
^^^^^
   jacquard translate <input_dir> <output_dir> [OPTIONS]

positional arguments:
  input                 Directory containing VCF files (and VarScan high
                        confidence files). Other file types ignored
  output                Directory containing VCF files. Will create if doesn't
                        exist and will overwrite files in output directory as
                        necessary

optional arguments:
  --allow_inconsistent_sample_sets
                        Set this flag if not every patient is represented by
                        the same set of caller-VCFs.
  --varscan_hc_filter_file_regex VARSCAN_HC_FILTER_FILE_REGEX
                        Regex pattern that identifies optional VarScan
                        high-confidence filter files. The VCF, high-confidence
                        file pairs should share the same prefix. For example,
                        given patientA.snp.vcf, patientA.indel.vcf,
                        patientA.snp.fpfilter.pass, and
                        patientA.indel.fpfilter.pass, you could enable this
                        option as varscan_hc_filter_file_regex='.fpfilter.pass$'


   VarScan tips:
   To translate VarScan calls, Jacquard requires the VarScan VCF files (snp
   and/or indel). Jacquard can additionally accept VarScan somatic high-
   confidence files; the translate command supresses VCF variant records that
   are absent in the high-confidence files. Somatic high-confidence files
   should be placed alongside corresponding VarScan VCFs and must have the same
   prefix as their corresponding VCF file. The suffix should be specified using
   the command line argument.

   Note: Optional high-confidence files are not VCF files.


Merge
^^^^^
The merge command takes an input directory of VCFs and filters both the
variants and format tags based on your specifications.


By default, merge contains only Jacquard-translated format tags (JQ_.*) and
includes all variants with valid syntax at loci where at least one variant was
somatic. The resulting filtered files contain fewer rows, yet higher quality
data than the input files.


After filtering, the merge command combines all of the input VCFs into a single
merged VCF that includes all necessary information for continuing your analysis.


The resulting VCF files contain the distinct set of all coordinates (CHROM, POS,
REF, and ALT) and samples from the input files, provided they pass the filters.
Each coordinate from the input VCF files is added to the output file, which
increases the file length. Additionally, sample columns are merged for each
patient, adding sample specific information and leading to increased column and
file width.


Importantly, rather than giving caller-wise sample columns in the output VCf
file, merge emits patient-wise sample columns. For each patient, the merge
command joins the set of corresponding sample columns into a single column. The
grouping of sample-specific information for each patient helps to easily
analyze the data.

   <merge diagram>

Usage
^^^^^
usage: jacquard merge <input_dir> <output_file> [OPTIONS]

positional arguments:
  input                 Directory containing VCF files. Other file types
                        ignored
  output                VCF file

optional arguments:
  --include_format_tags     Comma-separated user-defined list of regular
                            expressions for format tags to be included in output
  --include_cells           valid: Only include valid variants
                            all: Include all variants
                            passed: Only include variants which passed their
                                    respective filter
                            somatic: Only include somatic variants
  --include_rows            at_least_one_somatic: Include all variants at loci
                                                  where at least one variant
                                                  was somatic
                            all_somatic: Include all variants at loci where all
                                         variants were somatic
                            at_least_one_passed: Include all variants at loci
                                                 where at least one variant
                                                 passed
                            all_passed: Include all variants at loci where all
                                        variants passed
                            all: Include all variants at loci


Summarize
^^^^^^^^^
The summarize command uses the Jacquard-specific tags to aggregate caller
information from the file, providing a summary-level view. The inclusion of
summary fields, such as averages, helps you to easily determine which are the
true variants.

The summarized format tags contain the prefix 'JQ_SUMMARY'.

   <summarize diagram>

Usage
^^^^^
usage: jacquard summarize <input_file> <output_file>

positional arguments:
  input                Jacquard-merged VCF file (or any VCF with Jacquard
                       tags; e.g. JQ_SOM_MT)
  output               VCF file


Expand
^^^^^^
The expand command converts a VCF file into a tab-delimited file in a tabular
format. This format is more suitable than a VCF for analysis and visualization
in R, Pandas, Excel, or another third party application.

The 'fixed' fields (i.e. CHROM, POS, ID, REF, ALT, QUAL, FILTER) are directly
copied from the input VCF file. Based on the metaheaders, each field in the
INFO column is expanded into a separate column named after its tag ID. Also,
based on the metaheaders, each FORMAT tag is expanded into a set of columns,
one for each sample, named as <format tag ID>|<sample column name>.

This command also emits a tab-delimited glossary file, created based on the
metaheaders in the input VCF file. FORMAT and INFO tag IDs are listed in the
glossary and are defined by their metaheader description.

   <expand diagram: figure 3 - excel sheet & schematic from the poster>

Usage
^^^^^
usage: jacquard expand <input_file> <output_file> [OPTIONS]

positional arguments:
  input                 VCF file. Other file types ignored
  output                TXT file
optional arguments:
  -s SELECTED_COLUMNS_FILE, --selected_columns_file SELECTED_COLUMNS_FILE
                        File containing an ordered list
                        of column names to be included in the output file;
                        column names can include regular expressions.

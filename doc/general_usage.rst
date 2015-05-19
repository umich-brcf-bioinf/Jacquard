General usage
=============

usage: jacquard <SUBCOMMAND> [ARGUMENTS] [OPTIONS]

for help on a specific command:
 'jacquard <SUBCOMMAND> -h' 

subcommands: {translate, merge, summarize, expand} 

translate   Accepts a directory of VCF files (and VarScan high
confidence files). Creates a new directory of VCF files, adding Jacquard-specific
FORMAT tags to each VCF record. 

merge       Accepts a directory of VCF files and returns a single merged VCF file.

summarize   Accepts a Jacquard-tagged VCF file and creates a new file, adding
summary fields. 

expand      Converts a VCF file into an tab-delimited file in a tabular format 
which would be more suitable for analysis and visualization in R, Pandas, Excel 
or other third party applications. 

See https://github.com/umich-brcf-bioinf/Jacquard for more info

Jacquard first writes output files to a temporary directory and only copies them 
upon successful completion of each subcommand.

Error, warning, and info messages are written to console and log file. 
Debug messages are written to the log file unless logger is initialized as verbose (in
which case debug is also echoed to console). 

<General data flow diagram>

Input file conventions
^^^^^^^^^^^^^^^^^^^^^^
Jacquard assumes that the first element of the filename is a patient identifier.

   patientA-113.mutect.vcf patientA-113.strelka.snv.vcf
   patientA-113.strelka.indel.vcf

This set of three files all have the same patient identifier (patientA-113)
representing the same tumor-normal pair; these files will be combined into a
single pair of tumor-normals in the merged VCF. See merge for more details.


To translate a specific VCF dialect, Jacquard determines the source variant
caller based on metaheaders. For this reason it is essential that you preserve
all metaheaders from the source VCF.


For a specific source VCF, Jacquard automatically determines the normal and
tumor samples based on the column header and metaheaders.



Commands
========

Common options
^^^^^^^^^^^^^^
  -h, --help            Show this help message and exit 
  --force               Overwrite contents of output directory 
  --log_file LOG_FILE   Path to log file destination. Defaults to current 
                        working directory if not specified. 
  -v, --verbose         Prints detailed log to the console


translate
^^^^^^^^^
The translate command accepts a directory of VCF files and creates a new
directory of "translated" VCF files where translated VCFs are simply extended
to include several Jacquard-specific format tags.

You can gather all input VCFs into a single directory and run translate once, or
partition VCFs into separate directories (for example, by variant caller) and
run translate once for each input directory. When partitioning into separate
input directories, all file names must be unique.

The translated format tags contain a caller specific prefix; example: 'JQ_SK'
for Strelka, 'JQ_VS' for VarScan and 'JQ_MT' for MuTect.

   <translate diagram>

Usage
^^^^^
   jacquard translate [OPTIONS] <input_dir> <output_dir> 

positional arguments:
  input                 Directory containing VCF files (and VarScan high
  confidence files). Other file types ignored 
  output                Directory containing VCF files. Will create if doesn't exist and will overwrite files
  in output directory as necessary

optional arguments:
  --allow_inconsistent_sample_sets
                        Set this flag if not every
                        patient is represented by the same set of caller-VCFs.
  --varscan_hc_filter_file_regex VARSCAN_HC_FILTER_FILE_REGEX
                        Regex pattern that identifies
                        optional VarScan high-confidence filter files. The VCF,
                        high-confidence file pairs should share the same
                        prefix. For example, given patientA.snp.vcf,
                        patientA.indel.vcf, patientA.snp.fpfilter.pass,
                        patientA.indel.fpfilter.pass, you could enable this
                        option as varscan_hc_filter_file_regex='.fpfilter.pass$'


   VarScan tips To translate VarScan calls, Jacquard requires the VarScan VCF
   files (snp and/or indel). Jacquard can additionally accept VarScan somatic
   high-confidence files; the translate command deprecates VCF variant records
   absent in the high-confidence file. Somatic high-confidence files should be
   placed alongside corresponding VarScan VCFs and must have the same prefix as
   their corresponding VCF file. The suffix should be specified using the
   command line argument.
   
   Note: Optional high-confidence files are not VCF files.

Merge
^^^^^
The merge command takes an input directory of VCFs and filters the
variants and format tags based on your specifications.

By default merge contains only Jacquard-translated format tags
(JQ_.*) and includes all variants with valid syntax at loci where at least one variant was somatic. The resulting filtered files contain fewer rows, but higher quality data than the input files. 

The merge command then combines all the input VCFs into one merged VCF that
includes all necessary information for continuing your analysis.

The resulting VCF files contain all the coordinates and samples from the input
files, provided they pass the filters. It merges the sample columns for each
patient and adds sample specific information, leading to increased column and
file width. Also, each of the rows (coordinates) from the input VCF files are
added to the output file, which increases the file length. 

Usage
^^^^^

usage: jacquard merge <input_dir> <output_file> [OPTIONS]

positional arguments:
  input                 Directory containing VCF files. Other file types
  ignored output                VCF file

optional arguments:

  --include_format_tags     Comma-separated user-defined list of regular
  expressions for format tags to be included in output --include_cells         
  valid: Only include valid variants
                            all: Include all
                            variants passed: Only include variants which passed
                            their respective filter somatic: Only include
                            somatic variants
  --include_rows            at_least_one_somatic: Include all variants at loci
  where at least one variant was somatic
                            all_somatic: Include all
                            variants at loci where all variants were somatic
                            at_least_one_passed: Include all variants at loci
                            where at least one variant passed all_passed:
                            Include all variants at loci where all variants
                            passed all: Include all variants at loci

   <merge diagram>

Summarize
^^^^^^^^^
The summarize command uses the Jacquard-defined tags to aggregate the
caller information from the file, providing a summary view. The inclusion of
summary fields, such as averages, helps you to easily determine the true
variants.

The summarized format tags contain the prefix 'JQ_SUMMARY'.

Usage
^^^^^
usage: jacquard summarize <input_file> <output_file>

positional arguments:
  input                Jacquard-merged VCF file (or any VCF with Jacquard
  tags; e.g. JQ_SOM_MT) output               VCF file
  
   <summarize diagram>
   
Expand
^^^^^^
   
The expand command converts a VCF file into an tab-delimited file in a
tabular format which would be more suitable for analysis and visualization in
R, Pandas, Excel or other third party applications. 

   <expand diagram: figure 3 - excel sheet & schematic from the poster>
   
Usage
^^^^^
Usage: jacquard expand <input_file> <output_file>

positional arguments:
  input                 VCF file. Other file types ignored output             
  TXT file

optional arguments:
  -s SELECTED_COLUMNS_FILE, --selected_columns_file SELECTED_COLUMNS_FILE
                        File containing an ordered list
                        of column names to be included in the output file;
                        column names can include regular expressions.

   - The fixed fields are copied from the VCF directly. 
   - Based on the metaheaders, each INFO field is expanded into a separate
     column named after the INFO tag ID.
   - Based on the metaheaders, each FORMAT tag is expanded into a set of
     columns, one for each sample named as <FORMAT tag ID>|<sample column name>.
   

Note
^^^^
The Jacquard-produced output VCF is a fully compliant VCF file that can be 
easily loaded into an external program, such as an annotation tool.
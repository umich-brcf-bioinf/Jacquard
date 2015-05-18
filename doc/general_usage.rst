General usage
=============
<General data flow diagram>

Input file conventions
^^^^^^^^^^^^^^^^^^^^^^
Jacquard assumes that the first element of the filename is a patient identifier.

   patientA-113.mutect.vcf
   patientA-113.strelka.snv.vcf
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
  -h, --help            show this help message and exit
  --force               Overwrite contents of output directory
  --log_file LOG_FILE   Path to log file destination. Defaults to current working directory if not specified.
  -v, --verbose


translate
^^^^^^^^^
The translate command accepts a directory of VCF files and 
creates a new directory of "translated" VCF files where 
translated VCFs are simply extended to include several Jacquard-specific
format tags.

You can gather all input VCFs into a single directory and run translate once, 
or partition VCFs into separate directories (for example, by variant caller) 
and run translate once for each input directory. When partitioning into separate
input directories, all file names must be unique.

   <translate diagram>
Usage
^^^^^
   jacquard translate [OPTIONS] <input_dir> <output_dir> 

positional arguments:
  input                 Directory containing VCF files (and VarScan high confidence files). Other file types ignored
  output                Directory containing VCF files. Will create if doesn't exist and will overwrite files in output directory as necessary

optional arguments:
  --allow_inconsistent_sample_sets
                        Set this flag if not every patient is represented by the same set of caller-VCFs.
  --varscan_hc_filter_file_regex VARSCAN_HC_FILTER_FILE_REGEX
                        Regex pattern that identifies optional VarScan high-confidence filter files.
                        The VCF, high-confidence file pairs should share the same prefix.
                        For example, given patientA.snp.vcf, patientA.indel.vcf, patientA.snp.fpfilter.pass, patientA.indel.fpfilter.pass,
                        you could enable this option as varscan_hc_filter_file_regex='.fpfilter.pass$'


   VarScan tips
   To translate VarScan calls, Jacquard requires the VarScan VCF files 
   (snp and/or indel). Jacquard can additionally accept VarScan somatic 
   high-confidence files; the translate command deprecates VCF variant records
   absent in the high-confidence file. Somatic high-confidence files should be 
   placed alongside corresponding VarScan VCFs and must have the same prefix as 
   their corresponding VCF file. The suffix should be specified using the 
   command line argument.
   
   Note: Optional high-confidence files are not VCF files.
   
merge
   <merge diagram>
   including format tags
   filtering
summarize
   <summarize diagram>
expand
The expand command converts a VCF file into an tab-delimited file in a tabular 
format which would be more suitable for analysis and visualization in R, Pandas, 
Excel or other third party applications. 

   <expand diagram: figure 3 - excel sheet & schematic from the poster>
   
Usage: jacquard expand <input_file> <output_file>

positional arguments:
  input                 VCF file. Other file types ignored
  output                TXT file

optional arguments:
  -c COLUMN_SPECIFICATION, --column_specification COLUMN_SPECIFICATION
                        Path to text file containing column regular expressions to be included in output file

   - The fixed fields are copied from the VCF directly. 
   - Based on the metaheaders, each INFO field is expanded into a separated column.
   - 
   
General usage
   use of tmp dir
   logging
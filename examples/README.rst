=================
Example VCF Files
=================
The vcfs directory contains example VCF files to be used as input for merge.
These files have been de-identified and represent 5 patients.

Usage
=====

To use Jacquard with the example data, run the following commands. Note that
examples/vcfs is the input directory for translate.

**translate**

``$jacquard translate examples/00-input_vcfs <output_dir> [options]``

**merge**

``$jacquard merge examples/01-translated <output_vcf_file> [options]``

**summarize**

``$jacquard summarize examples/02-merged.vcf <output_vcf_file> [options]``

**expand**

``$jacquard expand exmaples/03-summarized.vcf <output_tsv_file> [options]``

|
|
Refer to Jacquard/README.rst for more usage details
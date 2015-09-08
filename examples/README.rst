=================
Example VCF Files
=================
The vcfs directory contains example VCF files to be used as input for merge.
These files have been de-identified and represent 5 patients.

Usage
=====

To use Jacquard with the example data, run the following commands. Note that
examples/vcfs is the input directory for translate.

*translate*
``$jacquard translate examples/vcfs <output_dir>/translated [options]``

*merge*
``$jacquard merge <output_dir>/translated <output_dir>/merged.vcf [options]``

*summarize*
``$jacquard summarize <output_dir>/merged.vcf <output_dir>/summarized.vcf [options]``

*expand*
``$jacquard expand <output_dir>/summarized.vcf <output_dir>/expanded.txt [options]``



Refer to Jacquard/README.rst for more usage details
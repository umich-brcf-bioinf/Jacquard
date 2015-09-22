===========
Quick Start
===========
This is a simple tutorial on how to use the four Jacquard commands.


1. **Install Jacquard (see** :ref:`installing-jacquard` **).**

|

2. **Unzip the** :download:`examples.zip <../examples/examples.zip>` **file to**
   **your home directory (or other directory of your choice).**

   The examples directory contains sample input VCFs from five patients run
   with three variant callers. The input VCFs are based on a subset of actual
   variant calls from clinical data; the samples were de-identified and VCF
   positions have been randomized to prevent downstream identification. As a
   result of randomization, the sample VCF reference calls at a specific
   position don't always match the base calls from a reference sequence.

   Along with the inputs, the example dierctory contains output from each
   Jacquard command, as explained below.
   
|

3. **Create an output directory.**

|

4. **The** :ref:`translate <translate-command>` **command creates new VCFs,**
   **adding a controlled vocabulary of new FORMAT tags.**
   ::

      $ jacquard translate examples/00-input_vcfs/ <output_dir>

|

5. **The** :ref:`merge <merge-command>` **command integrates a directory of**
   **VCFs into a single VCF.**
   ::

      $ jacquard merge examples/01-translated/ <output_vcf_file>

|

6. **The** :ref:`summarize <summarize-command>` **command adds new INFO fields**
   **and FORMAT tags that combine variant data from the merged VCF.**
   ::

      $ jacquard summarize examples/02-merged.vcf <output_vcf_file>

|

7. **The** :ref:`expand <expand-command>` **command explodes a VCF file into**
   **a tab-delimited file.**
   ::

      $ jacquard expand examples/03-summarized.vcf <output_tsv_file>


Refer to :ref:`Overview <overview-page>` for more information on Jacquard.

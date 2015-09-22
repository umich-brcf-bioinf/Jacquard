Frequently Asked Questions
==========================

**Is Jacquard a variant caller?**
   Jacquard is not a variant caller. It accepts VCF output from variant callers
   and integrates them for simplified annotation and analysis.



**Can Jacquard annotate data?**
   No, Jacquard cannot annotate data; however the output from *translate*,
   *merge*, and *summarize* can be run through an annotation tool such as
   `SnpEff <http://snpeff.sourceforge.net/index.html>`_ or
   `Annovar <http://annovar.openbioinformatics.org/en/latest/>`_.



**Can I use Jacquard with any variant caller?**
   *Merge* and *expand* are able to process VCF files from any variant caller.
   *Translate* and *summarize*, however, must be run with VCF files from one or
   more of the supported variant callers. Currently, Jacquard supports
   MuTect, VarScan, and Strelka.



**I'd like to merge my VCFs, but my caller isn't supported by Jacquard.**
   Both *merge* and *expand* commands can be used to show all of the results
   from different callers without standardization of the input data. However,
   it is recommended that the input data be standardized whenever possible to
   directly compare data across callers.



**Does Jacquard work with germline callers?**
   The *translate* command is optimized to work with tumor-normal sample pairs.
   Germline VCFs can be used with *merge* and *expand* commands. Better support
   for germline and pedigree VCFs is coming soon.


Still Have Questions?
^^^^^^^^^^^^^^^^^^^^^
Email bfx-jacquard@umich.edu for support and questions.

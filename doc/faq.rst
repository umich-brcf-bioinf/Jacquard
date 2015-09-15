Frequently Asked Questions
==========================

**Is Jacquard a variant caller?**

   No, Jacquard is not a variant caller.



**Can Jacquard annotate data?**

   No, Jacquard cannot annotate data; however the output from translate, merge,
   and summarize can be run through an annotation tool such as SnpEff or
   Annovar.



**Can I use Jacquard with any Variant Caller?**

   Merge and expand are able to process VCF files from any variant caller.
   Translate and summarize, however, must be run with VCF files from one or
   more of the supported variant callers. Currently, Jacquard supports
   MuTect, VarScan, and Strelka.



**Can I use Jacquard to show all the results from different callers without
standardization of the input data?**

   Both merge and expand can be used (either individually or together) to
   show all of the results from different callers without standardization of
   the input data. However, it is recommended that the input data be
   standardized whenever possible to directly compare data across callers.

Still Have Questions?
^^^^^^^^^^^^^^^^^^^^^
Email bfx-jacquard@umich.edu for support and questions.

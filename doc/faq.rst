Frequently Asked Questions
==========================

**Can I use Jacquard with other Variant Callers?**

   Merge and Expand are able to process VCF files from any variant caller.
   Translate and Summarize, however, must be run with VCF files from one or
   more of the supported variant callers. Currently, Jacquard supports
   MuTect, VarScan, and Strelka.



**Can I use Jacquard to show all the results from different callers without
standardization of the input data?**

   Both Merge and Expand can be used (either individually or together) to
   show all of the results from different callers without standardization of
   the input data. However, it is recommended that the input data be
   standardized whenever possible to directly compare data across callers.
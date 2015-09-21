Changelog
=========

0.42 (X/X/XXXX)
---------------
- Added docs on readthedocs.
- Improved workflow documentation with example data
- Merge will now disambiguate tag collisions from multiple VCs
- Translate/summarize now support GT tags
- Extended precision to 4 decimal places to support analysis of gene-panels.

0.41 (5/7/2015)
---------------
 - Combined filter command with merge command
 - Extended expand to create simple metaheader glossary
 - Adjusted code to support Python >=2.7 or 3.x
 - Improved checks for consistent VCF file sets
 - Fixed bug in merge that caused error if any VCFs were unsorted
 - Fixed bug in summarize that caused error if variant was called by subset of callers 

0.31 (3/17/2015)
----------------
 - Downgraded VCF format from 4.2 to 4.1
 - Fixed a bug that omitted CALLERS_REPORTED_LIST summary tag
 - Simplified summary tags; removed dependency on numpy
 - Adjusted VarScan translation to accept a file pattern to identify high-confidence files 


0.3 (3/9/2015)
--------------
 - Replaced [normalize], [tag] commands with [translate]; relaxed constraints on incoming data.
 - Renamed [consensus] to [summarize]
 - More consistent behavior in [expand]
 - Significantly improved [merge] performance 
 - Added new summary tags:
   - CALLERS_REPORTED_COUNT
   - CALLERS_REPORTED_LIST
   - SAMPLES_REPORTED_COUNT
   - CALLERS_PASSED_COUNT
   - CALLERS_PASSED_LIST
   - SAMPLES_PASSED_COUNT
 - Fixed bug in how Strelka calculated AF on indels
 - Improved command validation and error handling
 - Added project/code documentation 
 - Removed dependencies on pandas
  
  
0.21 (10/2014)
--------------
 - Initial public release



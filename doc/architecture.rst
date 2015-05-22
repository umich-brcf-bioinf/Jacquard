Architectural Overview
======================
This overview is intended for contributers.

Coding Standards & Guidelines
-----------------------------
 - Supports Python 2.7 and 3.x
 - Uses pylint to format the code
 - Separates words both in variable and method names with "_"
 - Separate words in the class names with CAPS
 - Uses absolute imports
 - Utilizes nosetests to run unit- and functional-tests

|

Commands
--------
Each command transforms files or directories and is indirectly executed through
the jacquard module.

Note that only jacquard.py is aware of all of the possible commands.

|

.. figure:: images/translate_uml_sequence.jpg

   **UML Sequence Diagram :** *An example UML sequence diagram for Translate.
   Other commands follow a similar, yet unique sequence.*


Translate
^^^^^^^^^
Translate standardizes VCF files by adding Jacquard-specific FORMAT tags.

*Filter flags are added to anomalous VCF records*
   Translate initializes several private classes that label anomalous records
   as such.

*New Jacquard-specific FORMAT tags are added to VCF records*
   Translate calls out to the utils/variant_caller_factory.py, which then adds
   new FORMAT tag values for each relevant variant caller.


Merge
^^^^^
Merge filters data from VCF files and then merges the files together.

*VCF records are filtered*
   Merge initializes a private class to filter the VCF records.

*VCF records are merged*
   Merge joins the VCF records together within its own module.


Summarize
^^^^^^^^^
Summarize aggregates sample-specific data based on Jacquard-specific FORMAT
tags.

*Sample-specific data is aggregated across callers*
   Summarize calls out to modules in utils/ to transform Jacquard-specific
   FORMAT tags into summarized FORMAT tags.


Expand
^^^^^^
Expand converts a VCF into a tab-delimited file in a tabular format.

*VCF fields are expanded*
   Expand separates INFO and FORMAT tags into distinct columns within its
   own module.

|

Variant Caller Transforms
-------------------------
Within this package are modules which transform VcfRecords. Each module
typically has a collection of tag classes, where each tag holds the metaheader
and code to transform a single VcfRecord.

Note that:

* A caller could manipulate any aspect of a VcfRecord, but (by strong
  convention) typically only adds information rather than deleting it. For
  example a sample-format tag, info field, filter field could be added.

* The only module in Jacuard that knows about all of the variant callers is
  variant_caller_factory.

|

Utils
-----
This package contains modules with methods that are relevant to multiple
commands.

vcf
^^^
The vcf module contains multiple classes that handle input and output files,
i.e., VcfReader and VcfRecord.

|

Test Conventions
----------------
Both unit- and functional-tests are written for Jaquard using the nosetests
framework.

Note that:

* Every method must have at least one corresponding test method
* Every class must have a TestCase
* Every command must have a functional test

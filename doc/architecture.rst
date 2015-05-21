Architectural Overview
======================
This overview is intended for contributers.

Coding Standards & Guidelines
-----------------------------
 - Uses pylint to format the code.
 - Uses '_' to separate the words both in variable and method names.
 - Uses caps to separate the words in the class names.
 - Uses absolute imports.
 - Supports Python 2.7 and 3.x.

Commands
--------
Jacquard is a suite of Python command line tools. Each tool is contained in a
module of the same name and called in the jacquard module. Only jacquard.py
is aware of all of the possible commands. 

Translate
^^^^^^^^^
Translate standardizes VCF files by adding Jacquard-specific format tags.

There are two main functions that Translate implements:

**1. Filter flags are added to anomalous VCF records.**
   Translate initializes several private classes that label anomalous records
   as such.

**2. New Jacquard-specific format tags are added to VCF records.**
   Translate calls out to the variant_caller_factory.py, which then adds new 
   format tag values for each relevant variant caller. Only
   variant_caller_factory.py is aware of all of the possible variant callers.


Merge
^^^^^



Summarize
^^^^^^^^^

Expand
^^^^^^

Transforms and Tags
-------------------

Utils
-----

Test Conventions
----------------

from collections import OrderedDict
import os
import re
import jacquard.jacquard_utils as jacquard_utils
import varscan, mutect, strelka, unknown

class TagProcessor():
    def __init__(self, callers):
        self.callers = callers
        
    def add_tags(self, input_line):
        no_newline_line = input_line.rstrip("\n")
        original_vcf_fields = no_newline_line.split("\t")
        new_vcf_fields = original_vcf_fields[:8]
        alt = original_vcf_fields[4]
        filter = original_vcf_fields[6]
        info = original_vcf_fields[7]
        format = original_vcf_fields[8]
        samples = original_vcf_fields[9:]

        count = 0
        for sample in samples:
            format_dict = OrderedDict(zip(format.split(":"), sample.split(":")))
            for tag in self.tags:
                format_dict = tag.format(alt, filter, info, format_dict, count)
                
            if count < 1: ##only add format column once
                new_vcf_fields.append(":".join(format_dict.keys()))
            new_vcf_fields.append(":".join(format_dict.values()))
            count += 1

        return "\t".join(new_vcf_fields) + "\n"
    
def execute(args, execution_context):
    pass
## Make TagProcessor initialize all the callers. Each caller will have a list "tags" as a parameter - input all the initialized tags for each.
## TagProcessor mimics LineProcessor in tags.
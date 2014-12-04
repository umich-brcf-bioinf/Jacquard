from __future__ import print_function
__version__ = 0.21

from collections import OrderedDict
from operator import itemgetter
import os
from os import listdir
import sys
import logger

global caller_versions
#TODO cgates: These should be in the callers or the caller factory, but not here.
caller_versions = {"VarScan":"v2.3", "MuTect": "v1.1.4", "Strelka": "v2.0.15"}

global jq_somatic_tag
global jq_af_tag
global jq_dp_tag
jq_somatic_tag = "HC_SOM"
jq_af_tag = "AF"
jq_dp_tag = "DP"


class JQException(Exception):
    def __init__(self, msg, *args):
        error_msg = msg.format(*[str(i) for i in args])
        super(JQException, self).__init__(error_msg)
        
    """Base class for exceptions in this module."""
    pass

def validate_directories(input_dir, output_dir):    
    if not os.path.isdir(input_dir):
        logger.error("Specified input directory [{}] does not exist.",
                     input_dir)
        exit(1)
    try:
        listdir(input_dir)
    except:
        logger.error("Specified input directory [{}] cannot be read. "+
            "Check permissions and try again.",input_dir)
        exit(1)

    if not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir)
        except:
            logger.error("Output directory [{}] could not be created. "+
                "Check parameters and try again", output_dir)
            exit(1)

def write_output(writer, headers, actual_sorted_variants):
    for line in headers:
        writer.write(line)
    for line in actual_sorted_variants:
        writer.write(line)

def sort_headers(headers):
    meta_headers = []
    field_header = ""
    for header in headers:
        if header.startswith("##"):
            header = header.replace("\t", "")
            meta_headers.append(header)
        else:
            field_header = header

    meta_headers.append(field_header)

    return meta_headers

def sort_data(all_variants):
    new_variants = []
    for variant in all_variants:
        split_variant = variant.split("\t")
        new_line = change_pos_to_int(split_variant)
        new_variants.append(new_line)

    #sort by CHROM, POS, REF, ALT
    sorted_variants = sorted(new_variants, key=itemgetter(0,1,3,4))

    actual_sorted_variants = []
    for variant in sorted_variants:
        new_field = [str(field) for field in variant]
        if "chr" not in new_field[0]:
            new_field[0] = "chr" + new_field[0]
        actual_sorted_variants.append("\t".join(new_field))

    return actual_sorted_variants

def change_pos_to_int(split_line):
    new_line = []
    split_line[0] = split_line[0].strip("chr")
    for field in split_line:
        try:
            new_field = int(field)
            new_line.append(new_field)
        except:
            new_line.append(field)
    return new_line

def combine_format_values(format, sample):
    return OrderedDict(zip(format.split(":"), sample.strip("\n").split(":")))

def combine_format_dict(format, sample):
    return OrderedDict(zip(format, sample))


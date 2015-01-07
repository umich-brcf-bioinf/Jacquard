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


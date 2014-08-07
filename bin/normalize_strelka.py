#!/usr/bin/python2.7
from collections import defaultdict
import glob
import os
from os import listdir
import re

import jacquard_utils
import normalize_utils

def add_subparser(subparser):
    parser_normalize_vs = subparser.add_parser("normalize_strelka", help="Accepts a directory containing Strelka VCF snvs/indel results and creates a new directory of merged, sorted VCFs")
    parser_normalize_vs.add_argument("input_dir", help="Path to directory containing VCFs. Each sample must have exactly two files matching these patterns: <sample>.indel.vcf, <sample>.snvs.vcf")
    parser_normalize_vs.add_argument("output_dir", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    jacquard_utils.validate_directories(input_dir, output_dir)

    normalize_utils.merge_and_sort(input_dir, output_dir, "strelka", execution_context)
#!/usr/bin/python2.7
from collections import defaultdict
import glob
import os
from os import listdir
import re

import jacquard_utils
import normalize_utils

def add_subparser(subparser):
    parser_normalize_vs = subparser.add_parser("normalize_varscan", help="Accepts a directory containing VarScan VCF snp/indel results and creates a new directory of merged, sorted VCFs with added high confidence tags")
    parser_normalize_vs.add_argument("input_dir", help="Path to directory containing VCFs. Each sample must have exactly two files matching these patterns: <sample>.indel.vcf, <sample>.snp.vcf, <sample>.indel.Germline.hc, <sample>.snp.Germline.hc, <sample>.indel.LOH.hc, <sample>.snp.LOH.hc, <sample>.indel.Somatic.hc, <sample>.snp.somatic.hc, <sample>.indel.*.hc, <sample>.snp.*.hc ")
    parser_normalize_vs.add_argument("output_dir", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    jacquard_utils.validate_directories(input_dir, output_dir)

    normalize_utils.merge_and_sort(input_dir, output_dir, "varscan", execution_context)
    
#!/usr/bin/python2.7
__version__ = 0.1

import argparse
import datetime
import numpy
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
from pandas import *
import random
import re
import sys 
import time

from pivot_variants import get_headers_and_readers, process_files, determine_input_keys
from tag_mutect import AlleleFreqTag, DepthTag, SomaticTag, LineProcessor, FileProcessor, tag_mutect_files
    
def validate_directories(input_dir, output_dir):    
    if not os.path.isdir(input_dir):
            print "Error. Specified input directory {0} does not exist".format(input_dir)
            exit(1)
            
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
            
if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pd.set_option('chained_assignment', None)

    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter, 
    description='''type 'Jacquard -h <subcommand>' for help on a specific subcommand''', 
    epilog="authors: Jessica Bene, Chris Gates 05/2014")
    parser.add_argument("-v", "--version", action='version', version="v{0}".format(__version__))
    subparsers = parser.add_subparsers(title="subcommands", dest="subparser_name")
    
    parser_pivot = subparsers.add_parser("pivot", help="Pivots input files so that given sample specific information is fielded out into separate columns. Returns an Excel file containing concatenation of all input files.")
    parser_pivot.add_argument("input_dir")
    parser_pivot.add_argument("output_file")
    parser_pivot.add_argument("-k", "--keys",
            help="Columns to be used as keys for the pivoting. Default keys for VCF are CHROM,POS,REF,ALT. Default keys for Epee TXT are CHROM,POS,REF,ANNOTATED_ALLELE,GENE_SYMBOL")
    parser_pivot.add_argument("-t", "--tags",
            help="Format tags to be fielded out in the pivoting.")
    
    parser_tagMutect = subparsers.add_parser("tagMutect", help="Accepts a directory of VCf results and creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser_tagMutect.add_argument("input_dir")
    parser_tagMutect.add_argument("output_dir")
    
    
    args = parser.parse_args()
    
    if args.subparser_name == "tagMutect":
        input_dir = os.path.abspath(args.input_dir)
        output_dir = os.path.abspath(args.output_dir)
        
        validate_directories(input_dir, output_dir)     
        input_metaheaders = ["##jacquard.version={0}".format(__version__), "##jacquard.tagMutect.command=tagMutect {0} {1}".format(input_dir, output_dir), "##jacquard.tagMutect.cwd={0}".format(os.path.dirname(os.getcwd()))]
        
        tag_mutect_files(input_dir, output_dir, input_metaheaders)
        
    elif args.subparser_name == "pivot":
        input_dir = os.path.abspath(args.input_dir)
        output_path = os.path.abspath(args.output_file)
        input_keys = args.keys.split(",") if args.keys else determine_input_keys(input_dir)
        pivot_values = args.tags.split(",") if args.tags else ["GT"]
        
        output_dir, outfile_name = os.path.split(output_path)
    
        validate_directories(input_dir, output_dir)
            
        fname, extension = os.path.splitext(outfile_name)
        if extension != ".xlsx": 
            print "Error. Specified output {0} must have a .xlsx extension".format(output_path)
            exit(1)
        
        sample_file_readers, headers, header_names, first_line = get_headers_and_readers(input_dir)
        process_files(sample_file_readers, input_dir, output_path, input_keys, pivot_values, headers, header_names, first_line)
    
    

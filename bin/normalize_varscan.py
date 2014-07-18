#!/usr/bin/python2.7
from collections import defaultdict
import glob
from operator import itemgetter, attrgetter
import os
from os import listdir
import re
# import vcf
# from vcf import utils

def identify_merge_candidates(in_files):
    merge_candidates = defaultdict(list)
    for in_file in in_files:
        merged_fname = re.sub("snp|indel", "merged", os.path.basename(in_file))
        merge_candidates[merged_fname].append(in_file)
        
    return merge_candidates

def get_headers(file):
    headers = []
    fh = open(file, "r")
    for line in fh:
        if line.startswith("#"):
            headers.append(line)
    fh.close()
    
    return headers

def validate_split_line(split_line, invalid):
    if re.search('\+|-|/', split_line[4]) or re.search('\+|-|/', split_line[3]):
        if not re.search('SS=5', split_line[7]):
            print "ERROR: {0}".format("\t".join(split_line).strip("\n"))
            invalid += 1
            
    return invalid
    
def change_pos_to_int(split_line):
    new_line = []
    for field in split_line:
        try:
            new_field = int(field)
            new_line.append(new_field)
        except:
            new_line.append(field)
    return new_line

def merge_data(files):
    all_variants = []
    invalid = 0
    for file in files:
        f = open(file, "r")
        for line in f:
            split_line = line.split("\t")

            if line.startswith("#"):
                continue
            else:
                invalid = validate_split_line(split_line, invalid)
                new_line = change_pos_to_int(split_line) #to sort properly
                all_variants.append(new_line)
        f.close()
        
    if invalid != 0:
        print "ERROR: {0} record(s) had illegal characters in REF or ALT. Review input files and try again".format(invalid)
        exit(1)
        
    return all_variants
                
def sort_data(all_variants):
    sorted_variants = sorted(all_variants, key=itemgetter(0,1,3,4)) #sort by CHROM, POS, REF, ALT
    
    variants = []
    for variant in sorted_variants:
        new_field = [str(field) for field in variant]
        variants.append("\t".join(new_field))
        
    return variants

def write_output(writer, headers, variants):
    for line in headers:
        writer.write(line)
    for line in variants:
        writer.write(line)

def merge(merge_candidates, output_dir):
    for merge_file in merge_candidates.keys():
        pair = merge_candidates[merge_file]
        headers = get_headers(pair[0])
        headers.append("##jacquard.normalize_varscan.sources={0},{1}\n".format(pair[0], pair[1]))
        
        all_variants = merge_data(pair)
        sorted_variants = sort_data(all_variants)
        
        merge_output = os.path.join(output_dir, merge_file)
        out_file = open(merge_output, "w")
        write_output(out_file, headers, sorted_variants)
        out_file.close()
        
    print "Wrote [{0}] VCF files to [{1}]".format(len(merge_candidates.keys()), output_dir)

def merge_and_sort(input_dir, output_dir, execution_context=[]):
    print "\n".join(execution_context)
    
    in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
    
    merge_candidates = identify_merge_candidates(in_files)
    
    total_files = 0
    for key, vals in merge_candidates.items():
        total_files += len(vals)
    print "Processing [{0}] samples from [{1}] files in [{2}]".format(len(merge_candidates.keys()), total_files, input_dir)
    
    merge(merge_candidates, output_dir)

def add_subparser(subparser):
    parser_normalize_vs = subparser.add_parser("normalize_varscan", help="Accepts a directory containing VarScan VCF snp/indel results and creates a new directory of merged, sorted VCFs")
    parser_normalize_vs.add_argument("input_dir", help="Path to directory containing VCFs. Each sample must have exactly two files matching these patterns: <sample>.indel.vcf, <sample>.snp.vcf")
    parser_normalize_vs.add_argument("output_dir", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")
    
def validate_directories(input_dir, output_dir):    
    print input_dir
    if not os.path.isdir(input_dir):
        print "Error. Specified input directory {0} does not exist".format(input_dir)
        exit(1)
    try:
        listdir(input_dir)
    except:
        print "Error: Specified input directory [{0}] cannot be read. Check permissions and try again.".format(input_dir)
        exit(1)
        
    if not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir)
        except:
            print "Error: Output directory could not be created. Check parameters and try again"
            exit(1)

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    validate_directories(input_dir, output_dir)

    merge_and_sort(input_dir, output_dir)
    
#!/usr/bin/python2.7
from collections import defaultdict
import glob
import os
from os import listdir
import re

import jacquard_utils

def identify_hc_variants(hc_candidates):
    hc_variants = {}
    for key, vals in hc_candidates.items():
        for file in vals:
            f = open(file, "r")
            for line in f:
                split_line = line.split("\t")
                if line.startswith("chrom"):
                    continue
                else:
                    hc_key = "^".join([split_line[0], split_line[1], split_line[2], split_line[3]])
                    hc_variants[hc_key] = 1
            
    return hc_variants

def write_to_merged_file(new_lines, headers, key):
    sorted_variants = jacquard_utils.sort_data(new_lines)
    
    writer = open(key, "w")
    jacquard_utils.write_output(writer, headers, sorted_variants)
    writer.close()

def mark_hc_variants(hc_variants, merge_candidates, output_dir):
    marked_as_hc = []

    for key, vals in merge_candidates.items():
        new_lines = []
        headers = []
        f = open(key, "r")
        for line in f:
            split_line = line.split("\t")
            if line.startswith('"'):
                line = line.replace("\t", "")
                line = line[1:-2] + "\n"
            if line.startswith("#"):
                headers.append(line)
            else:
                merge_key = "^".join([split_line[0], str(split_line[1]), split_line[3], split_line[4]])
                if merge_key in hc_variants:
                    if "JQ_HC_VS" not in split_line[7]:
                        split_line[7] += ";JQ_HC_VS"
                    marked_as_hc.append(merge_key)
                new_line = "\t".join(split_line)
                new_lines.append(new_line)
        f.close()
        
        sorted_headers = jacquard_utils.sort_headers(headers)
        write_to_merged_file(new_lines, sorted_headers, key)
    
    print "Wrote [{0}] VCF files to [{1}]".format(len(merge_candidates.keys()), output_dir)
    
    return marked_as_hc
    
def check_for_missing(required_vals, val, key, missing):
    missing_files = []
    for item in required_vals:
        if item not in val:
            missing_files.append(item)
            missing = 1
    if missing_files != []:
        print "ERROR: [{0}] is missing required files {1}".format(key.strip("_"), missing_files)
    
    return missing_files, missing

def check_for_unknown(required_vals, val, key, added):
    added_files = []
    for thing in val:
        if thing not in required_vals:
            added_files.append(thing)
            added = 1
    if added_files != []:
        print "WARNING: [{0}] has unknown .hc files {1}".format(key.strip("_"), added_files)
    
    return added_files, added
            
def validate_file_set(all_keys):
    sample_files = defaultdict(list)
    for key in all_keys:
        prefix = key.split("merged")[0]
        suffix = key.split("merged")[1]
        sample_files[prefix.strip("_")].append(suffix)

    required_vals = [".Germline.hc", ".LOH.hc", ".Somatic.hc", ".vcf"]
    missing = 0
    added = 0
    for key, val in sample_files.items():
        missing_files, missing = check_for_missing(required_vals, val, key, missing)
        added_files, added = check_for_unknown(required_vals, val, key, added)
        
    if missing == 1:
        print "ERROR: Some required files were missing. Review input directory and try again"
        exit(1)
    if added == 1:
        print "WARNING: Some samples had unknown .hc files"
        
    return sample_files
    
def identify_merge_candidates(in_files, out_dir):
    merge_candidates = defaultdict(list)
    hc_candidates = defaultdict(list)
    
    for in_file in in_files:
        fname, extension = os.path.splitext(in_file)
        if extension == ".vcf":
            merged_fname = re.sub("snp|indel", "merged", os.path.join(out_dir, os.path.basename(in_file)))
            merge_candidates[merged_fname].append(in_file)
        elif extension == ".hc":
            merged_fname = re.sub("snp|indel", "merged", os.path.join(out_dir, os.path.basename(in_file)))
            hc_candidates[merged_fname].append(in_file)
    
    all_keys = hc_candidates.keys() + merge_candidates.keys()
    sample_files = validate_file_set(all_keys)

    return merge_candidates, hc_candidates

def get_headers(file):
    headers = []
    fh = open(file, "r")
    for line in fh:
        if line.startswith("#"):
            headers.append(line)
    fh.close()
    
    return headers

def validate_split_line(split_line, invalid, warn):
    valid_line = []
    if re.search('\+|-|/', split_line[4]) or re.search('\+|-|/', split_line[3]):
        if not re.search('SS=5', split_line[7]):
            print "ERROR: {0}".format("\t".join(split_line).strip("\n"))
            invalid += 1
        elif re.search('SS=5', split_line[7]):
            warn += 1
    else:
        valid_line = split_line
            
    return (invalid, warn, valid_line)

def merge_data(files):
    all_variants = []
    invalid = 0
    warn = 0
    for file in files:
        f = open(file, "r")
        for line in f:
            split_line = line.split("\t")

            if line.startswith("#"):
                continue
            else:
                invalid, warn, valid_line = validate_split_line(split_line, invalid, warn)
                if valid_line != []:
                    all_variants.append("\t".join(valid_line))
        f.close()
    
    if warn != 0:
        print "WARNING: {0} records with SS=5 (somatic status is unknown) had illegal characters in REF or ALT. These lines will be excluded from the output.".format(warn)
    if invalid != 0:
        print "ERROR: {0} record(s) had illegal characters in REF or ALT. Review input files and try again".format(invalid)
        exit(1)
        
    return all_variants

def merge(merge_candidates, output_dir):
    for merge_file in merge_candidates.keys():
        pair = merge_candidates[merge_file]
        headers = get_headers(pair[0])
        
        sample_sources = "##jacquard.normalize_varscan.sources={0},{1}\n".format(os.path.basename(pair[0]), os.path.basename(pair[1]))
        print os.path.basename(merge_file) + ": " + sample_sources
        
        headers.append(sample_sources)
        headers.append('##INFO=<ID=JQ_HC_VS,Number=1,Type=Flag,Description="Jaquard high-confidence somatic flag for VarScan. Based on intersection with filtered VarScan variants.">\n')
        all_variants = merge_data(pair)
        sorted_variants = jacquard_utils.sort_data(all_variants)
        
        out_file = open(merge_file, "w")
        jacquard_utils.write_output(out_file, headers, sorted_variants)
        out_file.close()

def merge_and_sort(input_dir, output_dir, execution_context=[]):
    print "\n".join(execution_context)

    in_files = sorted(glob.glob(os.path.join(input_dir,"*")))
    
    merge_candidates, hc_candidates = identify_merge_candidates(in_files, output_dir)

    total_files = 0
    for key, vals in merge_candidates.items():
        total_files += len(vals)
    print "Processing [{0}] samples from [{1}] files in [{2}]".format(len(merge_candidates.keys()), total_files, input_dir)
    
    merge(merge_candidates, output_dir)

    hc_variants = identify_hc_variants(hc_candidates)
    marked_as_hc = mark_hc_variants(hc_variants, merge_candidates, output_dir)

def add_subparser(subparser):
    parser_normalize_vs = subparser.add_parser("normalize_varscan", help="Accepts a directory containing VarScan VCF snp/indel results and creates a new directory of merged, sorted VCFs with added high confidence tags")
    parser_normalize_vs.add_argument("input_dir", help="Path to directory containing VCFs. Each sample must have exactly two files matching these patterns: <sample>.indel.vcf, <sample>.snp.vcf, <sample>.indel.Germline.hc, <sample>.snp.Germline.hc, <sample>.indel.LOH.hc, <sample>.snp.LOH.hc, <sample>.indel.Somatic.hc, <sample>.snp.somatic.hc, <sample>.indel.*.hc, <sample>.snp.*.hc ")
    parser_normalize_vs.add_argument("output_dir", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    jacquard_utils.validate_directories(input_dir, output_dir)

    merge_and_sort(input_dir, output_dir, execution_context)
    
from collections import defaultdict
import glob
import os
import re
import shutil

import variant_callers.varscan 
import variant_callers.strelka
import variant_callers.mutect
import variant_callers.variant_caller_factory as variant_caller_factory

import utils as utils
import vcf as vcf

def get_headers(in_file):
    meta_headers = []
    header = ""
    fh = open(in_file, "r")
    for line in fh:
        if line.startswith("##"):
            meta_headers.append(line)
        elif line.startswith("#"):
            header = line
        else:
            break
    fh.close()
    
    return meta_headers, header

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
        
    return all_variants

def merge(merge_candidates, output_dir, caller):
    for merge_file in merge_candidates.keys():
        pair = merge_candidates[merge_file]
        meta_headers, header = get_headers(pair[0])
        sample_sources = caller.meta_header.format(os.path.basename(pair[0]), os.path.basename(pair[1]))
        print os.path.basename(merge_file) + ": " + sample_sources
        
        meta_headers.append(sample_sources)
        meta_headers.append(header)
        all_variants = merge_data(pair)

        sorted_variants = utils.sort_data(all_variants)
        
        out_file = open(merge_file, "w")
        utils.write_output(out_file, meta_headers, sorted_variants)
        out_file.close()

#TODO (cgates): Adjust per EX-91
def identify_merge_candidates(in_files, out_dir, caller):
    merge_candidates = defaultdict(list)
    hc_candidates = defaultdict(list)
    
    for in_file in in_files:
        #fname, extension = os.path.splitext(in_file)
        if in_file.lower().endswith(".vcf"):
            merged_fname = re.sub(caller.file_name_search, "merged", os.path.join(out_dir, os.path.basename(in_file)))
            merge_candidates[merged_fname].append(in_file)
        elif in_file.lower().endswith(".somatic.hc"):
            hc_candidates = caller.handle_hc_files(in_file, out_dir, hc_candidates)
    
    all_keys = hc_candidates.keys() + merge_candidates.keys()
    #caller.validate_file_set(all_keys)

    return merge_candidates, hc_candidates

def _validate_single_caller(filepaths, get_caller):
    callers = set()
    try:
        for filepath in filepaths:
            if filepath.endswith(".vcf"):
                reader = vcf.VcfReader(filepath, get_caller)
                caller = get_caller(reader.metaheaders,
                                    reader.column_header,
                                    filepath)
                callers.add(caller)
    except utils.JQException:
        raise utils.JQException("Problem normalizing the input VCF files. "
                                "Review input and try again.")

    if len(callers) > 1:
        raise utils.JQException("The input directory contains VCFs from "
                                "different variant callers. Review input and "
                                "try again.")
    else:
        return iter(callers).next()

def merge_and_sort(input_dir, output_dir, callers, execution_context=[]):
    print "\n".join(execution_context)

    in_files = sorted(glob.glob(os.path.join(input_dir,"*")))
    caller = _validate_single_caller(in_files, variant_caller_factory.get_caller)

    if caller.name == "MuTect":
        count = 0
        for in_file_name in in_files:
            fname, extension = os.path.splitext(in_file_name)
            if extension == ".vcf":
                shutil.copy(in_file_name, output_dir)
                count += 1
        print "Copied [{0}] VCF files from [{1}] to [{2}]".format(count, input_dir, output_dir)
        
    else:
        merge_candidates, hc_candidates = identify_merge_candidates(in_files, output_dir, caller)
    
        total_files = 0
        for key, vals in merge_candidates.items():
            total_files += len(vals)
        print "Processing [{0}] samples from [{1}] files in [{2}]".format(len(merge_candidates.keys()), total_files, input_dir)
        
        merge(merge_candidates, output_dir, caller)
        
        caller.final_steps(hc_candidates, merge_candidates, output_dir)
        
def add_subparser(subparser):
    parser_normalize_vs = subparser.add_parser("normalize", help="Accepts a directory containing VarScan VCF snp/indel results or Strelka VCF snvs/indels results and creates a new directory of merged, sorted VCFs with added high confidence tags")
    parser_normalize_vs.add_argument("input_dir", help="Path to directory containing VCFs. Each sample must have exactly two files matching these patterns: <sample>.indel.vcf, <sample>.snp.vcf, <sample>.indel.Germline.hc, <sample>.snp.Germline.hc, <sample>.indel.LOH.hc, <sample>.snp.LOH.hc, <sample>.indel.Somatic.hc, <sample>.snp.somatic.hc, <sample>.indel.*.hc, <sample>.snp.*.hc ")
    parser_normalize_vs.add_argument("output_dir", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    utils.validate_directories(input_dir, output_dir)

    callers = [variant_callers.strelka.Strelka(), variant_callers.varscan.Varscan(), variant_callers.mutect.Mutect()]
    merge_and_sort(input_dir, output_dir, callers, execution_context)
    
        

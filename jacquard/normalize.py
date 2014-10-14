from collections import defaultdict
import glob
import os
import re
import shutil
from utils import JQException, log
import collections
import vcf

from variant_callers import variant_caller_factory
import variant_callers.varscan 
import variant_callers.strelka
import variant_callers.mutect

import utils as utils


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

def identify_merge_candidates(in_files, out_dir, caller):
    merge_candidates = defaultdict(list)
    hc_candidates = defaultdict(list)
    
    for in_file in in_files:
        fname, extension = os.path.splitext(in_file)
        if extension == ".vcf":
            merged_fname = re.sub(caller.file_name_search, "merged", os.path.join(out_dir, os.path.basename(in_file)))
            merge_candidates[merged_fname].append(in_file)
        elif extension == ".hc":
            hc_candidates = caller.handle_hc_files(in_file, out_dir, hc_candidates)
    
    all_keys = hc_candidates.keys() + merge_candidates.keys()
    caller.validate_file_set(all_keys)

    return merge_candidates, hc_candidates

def validate_callers(in_file_name, validated_callers, callers):
    for caller in callers:
        in_file = open(in_file_name, "r")
        caller_name, valid = caller.validate_input_file(in_file)
        if valid == 1:
            if caller_name == "Unknown":
                print "ERROR: Unknown input file type"
                exit(1)
            if caller_name not in validated_callers:
                validated_callers.append(caller_name)
            break
        in_file.close()
    
    return validated_callers

def merge_and_sort(input_dir, output_dir, callers, execution_context=[]):
    print "\n".join(execution_context)

    in_files = sorted(glob.glob(os.path.join(input_dir,"*")))
    validated_callers = []
    for in_file_name in in_files:
        fname, extension = os.path.splitext(in_file_name)
        if extension == ".vcf":
            validated_callers = validate_callers(in_file_name, validated_callers, callers)

    if len(validated_callers) == 1:
        for initial_caller in callers:
            if validated_callers[0] == initial_caller.input_filepath:
                caller = initial_caller
    else:
        print "ERROR: It appears that input directory contains files called with different variant callers. Input directory should only have files for one caller. Review input and try again."
        exit(1)

    if caller.input_filepath == "MuTect":
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

def _partition_files_by_patient(input_dir, output_dir):
    in_files = sorted(glob.glob(os.path.join(input_dir, "*")))
    patient_to_files = defaultdict(list)
    for file_path in in_files:
        filename = os.path.basename(file_path)
        patient = filename.split(".")[0]
        patient_to_files[patient].append(vcf.FileReader(file_path))
    
    writer_to_readers = {}
    for patient in patient_to_files:
        file_writer = vcf.FileWriter(os.path.join(output_dir,patient+".normalized.vcf"))
        writer_to_readers[file_writer] = patient_to_files[patient]
        
    return writer_to_readers

def _determine_caller_per_directory(writer_to_readers, 
                                    get_caller=variant_caller_factory.get_caller):
    file_reader = writer_to_readers.values()[0][0]
    vcf_reader = vcf.VcfReader(file_reader)
    return get_caller(vcf_reader.metaheaders, vcf_reader.column_header, vcf_reader.file_name)

def _log_caller_info(vcf_readers):
    caller_count = collections.defaultdict(int)
    for vcf in vcf_readers:
        caller_count[vcf.caller.name] += 1
    for caller_name in sorted(caller_count):
        log("INFO: Recognized [{0}] {1} file(s)",
             caller_count[caller_name], caller_name)

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    utils.validate_directories(input_dir, output_dir)
    patient_to_files = _partition_files_by_patient(input_dir, output_dir)
    if not patient_to_files:
        log("ERROR: Specified input directory [{0}] contains no VCF files."
             "Check parameters and try again.", input_dir)
        #TODO cgates: move to jacquard.py
        shutil.rmtree(output_dir)
        exit(1)

    
    
    callers = [variant_callers.strelka.Strelka(), variant_callers.varscan.Varscan(), variant_callers.mutect.Mutect()]
    merge_and_sort(input_dir, output_dir, callers, execution_context)
    
        
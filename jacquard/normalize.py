from collections import defaultdict
import glob
import os
import re
import shutil
from utils import JQException
import collections

import variant_callers.varscan 
import variant_callers.strelka
import variant_callers.mutect
import variant_callers.variant_caller_factory as variant_caller_factory

import utils as utils
import vcf as vcf
import logger as logger

#TODO: hook up to something
def _validate_single_caller(filepaths, get_caller):
    callers = set()
    try:
        for filepath in filepaths:
            if filepath.endswith(".vcf"):
                reader = vcf.VcfReader(vcf.FileReader(filepath))
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

def add_subparser(subparser):
    parser = subparser.add_parser("normalize", help="Accepts a directory containing VarScan VCF snp/indel results or Strelka VCF snvs/indels results and creates a new directory of merged, sorted VCFs with added high confidence tags")
    parser.add_argument("input", help="Path to directory containing VCFs. Each sample must have exactly two files matching these patterns: <sample>.indel.vcf, <sample>.snp.vcf, <sample>.indel.Germline.hc, <sample>.snp.Germline.hc, <sample>.indel.LOH.hc, <sample>.snp.LOH.hc, <sample>.indel.Somatic.hc, <sample>.snp.somatic.hc, <sample>.indel.*.hc, <sample>.snp.*.hc ")
    parser.add_argument("output", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")
    
def _partition_input_files(in_files, output_dir, caller):
    patient_to_files = defaultdict(list)
        
    for file_path in in_files:
        basename = os.path.basename(file_path)
        patient = basename.split(".")[0]
        patient_to_files[patient].append(file_path)
   
    writer_to_readers = {}
    for patient,in_files in patient_to_files.items():
        output_file = caller.decorate_files(in_files, "normalized")
        file_writer = vcf.FileWriter(os.path.join(output_dir,output_file))
        file_readers = []
        
        for file_path in patient_to_files[patient]:
            file_readers.append(vcf.FileReader(file_path))
        
        writer_to_readers[file_writer] = file_readers

    return writer_to_readers

def _determine_caller_per_directory(in_files, 
                                    get_caller=variant_caller_factory.get_caller):
    for in_file in in_files:
        extension = os.path.splitext(os.path.basename(in_file))[1]
        if extension == ".vcf":
            file_reader = vcf.FileReader(in_file)
            vcf_reader = vcf.VcfReader(file_reader)

            return get_caller(vcf_reader.metaheaders, vcf_reader.column_header, vcf_reader.file_name)

def _log_caller_info(vcf_readers):
    caller_count = collections.defaultdict(int)
    
    for vcf in vcf_readers:
        caller_count[vcf.caller.name] += 1
    
    for caller_name in sorted(caller_count):
        logger.info("Recognized [{0}] {1} file(s)",
             caller_count[caller_name], caller_name)
        
def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    
    utils.validate_directories(input_dir, output_dir)
    in_files = sorted(glob.glob(os.path.join(input_dir, "*")))
    
    caller = _determine_caller_per_directory(in_files)
    logger.info("Recognized caller as {}", caller.name)
    
    caller.validate_vcfs_in_directory(in_files)    
    
    writer_to_readers = _partition_input_files(in_files, output_dir, caller)
    
    if not writer_to_readers:
        logger.error("Specified input directory [{0}] contains no files."
             "Check parameters and try again.", input_dir)
        
        #TODO cgates: move to jacquard.py
        shutil.rmtree(output_dir)
        exit(1)
        
    logger.info("Normalizing input files")
    count = 1
    for writer, readers in writer_to_readers.items():
        logger.info("Writing to [{}] ({}/{})", writer.output_filepath, count,
                    len(writer_to_readers))

        caller.normalize(writer, readers)
        count += 1
     

from __future__ import print_function
import glob
import os
#import shutil

from variant_callers import variant_caller_factory
import jacquard_utils
from jacquard_utils import JQException


def log(msg, *args):
    print(msg.format([str(i) for i in args]))

def write_headers(reader, writer, in_file, execution_context, caller):
    metaheader = ''.join(["{}\n".format(header) for header in execution_context])
    metaheader = caller.update_metaheader(metaheader)
    writer.write(in_file + "\n")
    for line in reader:
        if line.startswith("##"):
            writer.write(line)
        elif line.startswith("#"):
            writer.write(metaheader)
            writer.write("##jacquard.tag.caller={0}\n".format(caller))
            writer.write(line)
        else:
            writer.write(caller.add_tags(line))    
        
def add_subparser(subparser):
    parser_tag = subparser.add_parser("tag", help="Accepts a directory of VCf results and creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser_tag.add_argument("input_dir", help="Path to directory containing VCFs. Other file types ignored")
    parser_tag.add_argument("output_dir", help="Path to Jacquard-tagged VCFs. Will create if doesn't exist and will overwrite files in output directory as necessary")

def _get_header(input_dir, in_file):
    header = []
    with open(os.path.join(input_dir, in_file), "r") as vcf:
        for line in vcf:
            if line.startswith("##"):
                header.append(line)
            else:
                break
    return header
        

def determine_callers2(vcfProviders, factory):
    file_to_caller = {}
    any_file_failed = False
    for vcf in vcfProviders:
        try:    
            file_to_caller[vcf.name] = factory.get_caller(vcf)
        except factory.JQException as e:
            log("ERROR: Problem parsing [{}]:{}",vcf.name, e)
            any_file_failed = True
    if any_file_failed:
        raise JQException("There were problems parsing some VCFs. Review input files above and try again.")

    return file_to_caller



def determine_callers(input_dir, in_files, factory):
    file_to_caller = {}   
    problem_files = {"unknown":[],"invalid":[]}         
    for in_file in in_files:
        header = _get_header(input_dir, in_file)
        try:    
            file_to_caller[in_file] = factory.get_caller(header)
        except factory.JQCallerNotRecognizedException:
            problem_files["unknown"].append(in_file)
        except JQException:
            problem_files["invalid"].append(in_file)
    if len(problem_files["unknown"])>0:
        raise JQException("Unable to determine which caller was used on {0}. Check input files and try again.".format(problem_files["unknown"]))
    if len(problem_files["invalid"])>0:
        raise JQException("Invalid caller on {0}. Check input files and try again.".format(problem_files["invalid"]))
    return file_to_caller
    
def print_callers(file_to_caller):
    for key, vals in file_to_caller.items():
        print "Recognized [{0}] {1} file(s)".format(len(vals), key)
    
def tag_files(input_dir, output_dir, execution_context=[]):
    factory = variant_caller_factory.Factory()
    in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
    if len(in_files) < 1:
        print "ERROR. Specified input directory [{0}] contains no VCF files. Check parameters and try again."
        shutil.rmtree(output_dir)
        exit(1)

    print "\n".join(execution_context)
    print "Processing [{0}] VCF file(s) from [{1}]".format(len(in_files), input_dir)
    
    file_to_caller = determine_callers(input_dir, in_files, factory)
    print_callers(file_to_caller)

#     processors = {"VarScan" : FileProcessor(output_dir, tags=[varscan.AlleleFreqTag(), varscan.DepthTag(), varscan.SomaticTag()], execution_context_metadataheaders=execution_context), "MuTect": FileProcessor(output_dir, tags=[mutect.AlleleFreqTag(), mutect.DepthTag(), mutect.SomaticTag()], execution_context_metadataheaders=execution_context), "Strelka": FileProcessor(output_dir, tags=[strelka.AlleleFreqTag(), strelka.DepthTag(), strelka.SomaticTag()], execution_context_metadataheaders=execution_context)}
    
    total_number_of_files = len(file_to_caller)
    for count,filename,caller in enumerate(file_to_caller.items()):
        print "Reading [{0}] ({1}/{2})".format(os.path.basename(filename), count, total_number_of_files)
        fname, extension = os.path.splitext(os.path.basename(filename))
        new_file = fname + ".jacquardTags" + extension
        
        reader = open(os.path.join(input_dir, filename), "r")
        writer = open(os.path.join(output_dir, new_file), "w")
        write_headers(reader, writer, output_dir, os.path.basename(filename), execution_context, caller)
        
        reader.close()
        writer.close()
    
    out_files = sorted(glob.glob(os.path.join(output_dir,"*.vcf")))
    
    print "Wrote [{0}] VCF file(s) to [{1}]".format(len(out_files), output_dir)

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    jacquard_utils.validate_directories(input_dir, output_dir)
    
    tag_files(input_dir, output_dir, execution_context)
    
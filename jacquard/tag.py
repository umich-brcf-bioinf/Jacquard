from collections import OrderedDict, defaultdict
import glob
import os
import shutil

from variant_callers import variant_caller_factory
import jacquard_utils
from jacquard_utils import JQException


# class LineProcessor():
#     def __init__(self, tags):
#         self.tags = tags
#         
#     def add_tags(self, input_line):
#         no_newline_line = input_line.rstrip("\n")
#         original_vcf_fields = no_newline_line.split("\t")
#         new_vcf_fields = original_vcf_fields[:8]
#         alt = original_vcf_fields[4]
#         filter = original_vcf_fields[6]
#         info = original_vcf_fields[7]
#         format = original_vcf_fields[8]
#         samples = original_vcf_fields[9:]
# 
#         count = 0
#         for sample in samples:
#             format_dict = OrderedDict(zip(format.split(":"), sample.split(":")))
#             for tag in self.tags:
#                 format_dict = tag.format(alt, filter, info, format_dict, count)
#                 
#             if count < 1: ##only add format column once
#                 new_vcf_fields.append(":".join(format_dict.keys()))
#             new_vcf_fields.append(":".join(format_dict.values()))
#             count += 1
# 
#         return "\t".join(new_vcf_fields) + "\n"
# 
# class FileProcessor():
#     def __init__(self, output_dir, tags=[], execution_context_metadataheaders = []):
#         self._tags = tags
#         self._metaheader = self._metaheader_handler(execution_context_metadataheaders)
#         self._output_dir = output_dir
#         
#         for tag in self._tags:
#             self._metaheader += tag.metaheader
#         self._lineProcessor = LineProcessor(self._tags)
#             
#     def _metaheader_handler(self, metaheaders):
#         new_headers = ["{}\n".format(header) for header in metaheaders]
#         return ''.join(new_headers)
# 
#     def process(self, reader, writer, in_file_name, caller, handlers):
#         writer.write(handlers[in_file_name] + "\n")
#         for line in reader:
#             if line.startswith("##"):
#                 writer.write(line)
#             elif line.startswith("#"):
#                 if caller == "VarScan":
#                     if "NORMAL" in line and "TUMOR" in line:
#                         writer.write(self._metaheader)
#                         writer.write("##jacquard.tag.caller={0}\n".format(caller))
#                         writer.write(line)
#                     else:
#                         print "Unexpected VarScan VCF structure - missing NORMAL\t and TUMOR\n headers."
#                         shutil.rmtree(self._output_dir)
#                         exit(1)
#                 else:
#                     writer.write(self._metaheader)
#                     writer.write("##jacquard.tag.caller={0}\n".format(caller))
#                     writer.write(line)
#             else:
#                 edited_line = self._lineProcessor.add_tags(line)
#                 writer.write(edited_line)

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
        
def determine_caller(input_dir, in_files, factory):
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
        #TODO: Continue with getting tags starting from here!
#             valid = caller.validate_input_file(in_file)
#             if valid == 1:
#                 file_types[caller.name].append(file)
#                 if (caller.good):
#                     handler = "{0}: ##jacquard.tag.handler={1}".format(os.path.basename(file), caller.name)
#                     print handler
#                     handlers[os.path.basename(file)] = "##jacquard.tag.handler={0}".format(caller.name)
#                     inferred_caller = "##jacquard.tag.caller={0}".format(caller.name)                 
#                     print "{0}: {1}".format(os.path.basename(file), inferred_caller)                    
#                 else:  
#                     print "ERROR. {0}: ##jacquard.tag.handler={1}".format(os.path.basename(file), caller.name)
#                 break
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
    
    file_to_caller = determine_caller(input_dir, in_files, factory)
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
    
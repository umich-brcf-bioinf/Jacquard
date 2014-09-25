from collections import OrderedDict, defaultdict
import glob
import os
import shutil

from variant_callers import variant_caller_factory
import jacquard_utils


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

def file_processor(reader, writer, output_dir, in_file, handlers, execution_context, caller):
    metaheader = ''.join(["{}\n".format(header) for header in execution_context])
    metaheader = caller.update_metaheader(metaheader)
    writer.write(handlers[in_file] + "\n")
    for line in reader:
        if line.startswith("##"):
            writer.write(line)
        elif line.startswith("#"):
            if caller.validate_record(line):
                writer.write(metaheader)
                writer.write("##jacquard.tag.caller={0}\n".format(caller))
                writer.write(line)
            else:
                shutil.rmtree(output_dir)
                exit(1)
        else:
            writer.write(caller.add_tags(line))    
        
def add_subparser(subparser):
    parser_tag = subparser.add_parser("tag", help="Accepts a directory of VCf results and creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser_tag.add_argument("input_dir", help="Path to directory containing VCFs. Other file types ignored")
    parser_tag.add_argument("output_dir", help="Path to Jacquard-tagged VCFs. Will create if doesn't exist and will overwrite files in output directory as necessary")

def determine_file_types(input_dir, in_files, factory):
    file_types = defaultdict(list)   
    handlers = {}                              
    for _ in in_files:
        for caller in factory.callers:
            in_file = open(os.path.join(input_dir, file), "r")
            valid = caller.validate_input_file(in_file)
            if valid == 1:
                file_types[caller.name].append(file)
                if (caller.good):
                    handler = "{0}: ##jacquard.tag.handler={1}".format(os.path.basename(file), caller.name)
                    print handler
                    handlers[os.path.basename(file)] = "##jacquard.tag.handler={0}".format(caller.name)
                    inferred_caller = "##jacquard.tag.caller={0}".format(caller.name)                 
                    print "{0}: {1}".format(os.path.basename(file), inferred_caller)                    
                else:  
                    print "ERROR. {0}: ##jacquard.tag.handler={1}".format(os.path.basename(file), caller.name)
                break
    return file_types, handlers
    
def print_file_types(output_dir, file_types):
    for key, vals in file_types.items():
        print "Recognized [{0}] {1} file(s)".format(len(vals), key)
    if "Unknown" in file_types.keys():
        print "ERROR. Unable to determine which caller was used on [{0}]. Check input files and try again.".format(file_types["Unknown"])
        shutil.rmtree(output_dir)
        exit(1)
    
def tag_files(input_dir, output_dir, execution_context=[]):
    factory = variant_caller_factory.Factory()
    in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
    if len(in_files) < 1:
        print "ERROR. Specified input directory [{0}] contains no VCF files. Check parameters and try again."
        shutil.rmtree(output_dir)
        exit(1)

    print "\n".join(execution_context)
    print "Processing [{0}] VCF file(s) from [{1}]".format(len(in_files), input_dir)
    
    file_types, handlers = determine_file_types(input_dir, in_files, factory)
    print_file_types(output_dir, file_types)

#     processors = {"VarScan" : FileProcessor(output_dir, tags=[varscan.AlleleFreqTag(), varscan.DepthTag(), varscan.SomaticTag()], execution_context_metadataheaders=execution_context), "MuTect": FileProcessor(output_dir, tags=[mutect.AlleleFreqTag(), mutect.DepthTag(), mutect.SomaticTag()], execution_context_metadataheaders=execution_context), "Strelka": FileProcessor(output_dir, tags=[strelka.AlleleFreqTag(), strelka.DepthTag(), strelka.SomaticTag()], execution_context_metadataheaders=execution_context)}
    
    total_number_of_files = len(in_files)
    for count,in_file in enumerate(in_files):
        print "Reading [{0}] ({1}/{2})".format(os.path.basename(in_file), count, total_number_of_files)
        fname, extension = os.path.splitext(os.path.basename(in_file))
        new_file = fname + ".jacquardTags" + extension
        
        in_file_reader = open(os.path.join(input_dir, in_file), "r")
        out_file_writer = open(os.path.join(output_dir, new_file), "w")
        
        for key, vals in file_types.items():
            if in_file in vals:
                for caller in factory.callers:
                    if caller == key:
                        file_processor(in_file_reader, out_file_writer, output_dir, os.path.basename(in_file), handlers, execution_context, caller)
        
        in_file_reader.close()
        out_file_writer.close()
    
    out_files = sorted(glob.glob(os.path.join(output_dir,"*.vcf")))
    
    print "Wrote [{0}] VCF file(s) to [{1}]".format(len(out_files), output_dir)

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    jacquard_utils.validate_directories(input_dir, output_dir)
    
    tag_files(input_dir, output_dir, execution_context)
    
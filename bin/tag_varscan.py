#!/usr/bin/python2.7
from collections import OrderedDict
import glob
import os
from os import listdir
from os.path import isfile, join
import shutil

class AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_AF_VS,Number=1,Type=Float, Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ).">\n'

    def format(self, info_string, format_param_string, format_value_string):
        format_param_array = format_param_string.split(":")
        format_value_array = format_value_string.split(":")
        format_dict = dict(zip(format_param_array, format_value_array))
        
        if "FREQ" in format_dict.keys():
            format_value_string += ":" + self.roundTwoDigits(format_dict["FREQ"])
            format_param_string += ":JQ_AF_VS"
            
        return format_param_string, format_value_string

    def roundTwoDigits(self, value): 
        if len(value.split(".")[1]) <= 2:
            return value
        else:
            return str(round(100 * float(value))/100) 
        
class DepthTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_DP_VS,Number=1,Type=Float, Description="Jacquard depth for VarScan (based on DP).">\n'

    def format(self, info_string, format_param_string, format_value_string):
        format_param_array = format_param_string.split(":")
        format_value_array = format_value_string.split(":")
        format_dict = dict(zip(format_param_array, format_value_array))
        
        if "DP" in format_dict.keys():
            format_value_string += ":" + format_dict["DP"]
            format_param_string += ":JQ_DP_VS"
            
        return format_param_string, format_value_string
    
class SomaticTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_SOM_VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1= somatic (based on SOMATIC info tag and if sample is TUMOR).">\n'
#  
    def format(self, info_string, format_param_string, format_value_string):
        info_array = info_string.split(";")
#         format_param_array = format_param_string.split(":")
#         format_value_array = format_value_string.split(":")
#         format_dict = dict(zip(format_param_array, format_value_array))

        if "SS=2" in info_array:
            format_value_string += ":1"
            format_param_string += ":JQ_SOM_VS"
        return format_param_string, format_value_string
#  
#     def somatic_status(self, ss_value):
#         if ss_value == "2":
#             return "1"
#         else:
#             return "0"

class LineProcessor():
    def __init__(self, tags):
        self.tags = tags
        
    def add_tags(self, input_line):
        no_newline_line = input_line.rstrip("\n")
        original_vcf_fields = no_newline_line.split("\t")
        new_vcf_fields = original_vcf_fields[:8]
        info = original_vcf_fields[7]
        format = original_vcf_fields[8]
        samples = original_vcf_fields[9:]

        count = 0
        for sample in samples:
            format_dict = OrderedDict()
            for tag in self.tags:
                param, value = tag.format(info, format, sample)
                param_list = param.split(":")
                value_list = value.split(":")
                for i in range(0, len(param_list)):
                    format_dict[param_list[i]] = value_list[i]

            if count < 1: ##only add format column once
                new_vcf_fields.append(":".join(format_dict.keys()))
            new_vcf_fields.append(":".join(format_dict.values()))
            count += 1
        print new_vcf_fields 
        return "\t".join(new_vcf_fields) + "\n"

class FileProcessor():
    
    def __init__(self, tags=[], execution_context_metadataheaders = []):
        self._tags = tags
        self._metaheader = self._metaheader_handler(execution_context_metadataheaders)
        
        for tag in self._tags:
            self._metaheader += tag.metaheader
        self._lineProcessor = LineProcessor(self._tags)
            
    def _metaheader_handler(self, metaheaders):
        new_headers = ["##{}\n".format(header) for header in metaheaders]
        return ''.join(new_headers)

    def process(self, reader, writer):

        for line in reader:
            if line.startswith("##"):
                writer.write(line)
            elif line.startswith("#"):
                writer.write(self._metaheader)
                writer.write(line)
            else:
                edited_line = self._lineProcessor.add_tags(line)
                writer.write(edited_line)


def add_subparser(subparser):
    parser_tagVarscan = subparser.add_parser("tag_varscan", help="Accepts a directory of VCf results and creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser_tagVarscan.add_argument("input_dir")
    parser_tagVarscan.add_argument("output_dir")

def validate_directories(input_dir, output_dir):    
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
        
def tag_varscan_files(input_dir, output_dir, input_metaheaders=[]):
    processor = FileProcessor(tags=[AlleleFreqTag(), DepthTag(), SomaticTag()], execution_context_metadataheaders=input_metaheaders)
    
    in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
    if len(in_files) < 1:
        print "Error: Specified input directory [{0}] contains no VCF files. Check parameters and try again."
        exit(1)
        
    print "\n".join(input_metaheaders)
    print "Processing [{0}] VCF file(s) from [{1}]".format(len(in_files), input_dir)
                                                       
    for file in in_files:
        fname, extension = os.path.splitext(os.path.basename(file))
        if not isfile(file):
            continue
        new_file = fname + "_jacquard" + extension
        
        in_file = open(os.path.join(input_dir, file), "r")
        
        for line in in_file:
            if line.startswith("##source=VarScan2"):
                in_file.close()
                in_file = open(os.path.join(input_dir, file), "r")
                out_file = open(os.path.join(output_dir, new_file), "w")
                processor.process(in_file, out_file)
                out_file.close()
                break
        in_file.close()
    out_files = sorted(glob.glob(os.path.join(output_dir,"*.vcf")))
    
    print "Wrote [{0}] VCF file(s) to [{1}]".format(len(out_files), output_dir)

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
     
    validate_directories(input_dir, output_dir)
    tag_varscan_files(input_dir, output_dir, execution_context)
    
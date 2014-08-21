#!/usr/bin/python2.7
from collections import OrderedDict, defaultdict
import glob
import os
from os import listdir
from os.path import isfile, join
import shutil
import sys

import jacquard_utils

class Varscan():
    def __init__(self):
        self.name = "VarScan"
        
    def validate_input_file(self, input_file):
        valid = 0
        for line in input_file:
            if line.startswith("##source=VarScan2"):
                valid = 1
            elif line.startswith("##"):
                continue
            else:
                break
        return (self.name, valid)
        
class Mutect():
    def __init__(self):
        self.name = "MuTect"
        
    def validate_input_file(self, input_file):
        valid = 0
        for line in input_file:
            if line.startswith("##MuTect"):
                valid = 1
            elif line.startswith("##"):
                continue
            else:
                break
        return (self.name, valid)
     
class Strelka():
    def __init__(self):
        self.name = "Strelka"
        
    def validate_input_file(self, input_file):
        valid = 0
        for line in input_file:
            if line.startswith("##source=strelka"):
                valid = 1
            elif line.startswith("##"):
                continue
            else:
                break
        return (self.name, valid)
     
class Unknown():
    def __init__(self):
        self.name = "Unknown"
        
    def validate_input_file(self, input_file):
        valid = 1
        return (self.name, valid)
    
    
class Varscan_AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_AF_VS,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)

    def format(self, alt, filter, info_string, format_dict, count):
        if "FREQ" in format_dict.keys():
            freq = format_dict["FREQ"].split(",")
            format_dict["JQ_AF_VS"] = self.roundTwoDigits(freq)
            
        return format_dict

    def roundTwoDigits(self, value): 
        new_values = []
        for val in value:
            new_val = str(float(val.strip("%"))/100)
            if len(new_val.split(".")[1]) <= 2:
                new_values.append(new_val)
            else:
                new_values.append(str(round(100 * float(new_val))/100))
        return ",".join(new_values) 
        
class Varscan_DepthTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_DP_VS,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)

    def format(self, alt, filter, info_string, format_dict, count):
        if "DP" in format_dict.keys():
            format_dict["JQ_DP_VS"] = format_dict["DP"]

        return format_dict
    
class Varscan_SomaticTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_SOM_VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1= somatic (based on SOMATIC info tag and if sample is TUMOR),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)
#  
    def format(self, alt, filter, info_string, format_dict, count):
        info_array = info_string.split(";")

        if "SS=2" in info_array and "JQ_HC_VS" in info_array:
            format_dict["JQ_HC_SOM_VS"] = self.somatic_status(count)
        else:
            format_dict["JQ_HC_SOM_VS"] = "0"
            
        return format_dict
#  
    def somatic_status(self, count):
        if count == 0: #it's NORMAL
            return "0"
        else: #it's TUMOR
            return "1"

class Mutect_AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_AF_MT,Number=A,Type=Float,Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)

    def format(self, alt, filter, info, format_dict, count):
        if "FA" in format_dict.keys():
            freq = format_dict["FA"].split(",")
            format_dict["JQ_AF_MT"] = self.roundTwoDigits(freq)
            
        return format_dict

    def roundTwoDigits(self, value): 
        new_values = []
        for val in value:
            if len(val.split(".")[1]) <= 2:
                new_values.append(val)
            else:
                new_values.append(str(round(100 * float(val))/100))
        return ",".join(new_values)
        
class Mutect_DepthTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_DP_MT,Number=1,Type=Float,Description="Jacquard depth for MuTect (based on DP),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)

    def format(self, alt, filter, info, format_dict, count):
        if "DP" in format_dict.keys():
            format_dict["JQ_DP_MT"] = format_dict["DP"]
            
        return format_dict
    
class Mutect_SomaticTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_SOM_MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)

    def format(self, alt, filter, info, format_dict, count):
        if "SS" in format_dict.keys():
            format_dict["JQ_HC_SOM_MT"] = self.somatic_status(format_dict["SS"])
        else:
            format_dict["JQ_HC_SOM_MT"] = "0"
        return format_dict

    def somatic_status(self, ss_value):
        if ss_value == "2":
            return "1"
        else:
            return "0"

class Strelka_AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_AF_SK,Number=A,Type=Float,Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)

    def format(self, alt, filter, info, format_dict, count):
        afs = []
        if alt == ".":
            afs = ["."]
        else:
            split_alt = alt.split(",")
            for alt_allele in split_alt:
                if "AU" in format_dict.keys(): #if it's an snv
                    numerator = float(format_dict[alt_allele + "U"].split(",")[1])
                    tags = ["AU", "CU", "TU", "GU"]
                    denominator = 0
                    for tag in tags:
                        denominator += float(format_dict[tag].split(",")[1])
                    af = numerator/denominator if denominator != 0 else 0.0
                    
                elif "TAR" in format_dict.keys(): #if it's an indel
                    numerator = float(format_dict["TAR"].split(",")[1])
                    denominator = float(format_dict["DP2"])
                    af = numerator/denominator if denominator != 0 else 0.0
                else:
                    continue
                
                rounded_af = self.roundTwoDigits(str(af))
                capped_af = min(rounded_af, "1.00")
                afs.append(capped_af)
        
        if afs != []:
            format_dict["JQ_AF_SK"] = ",".join(afs)
           
        return format_dict

    def roundTwoDigits(self, value): 
        if len(value.split(".")[1]) <= 2:
            return value
        else:
            return str(round(100 * float(value))/100) 

class Strelka_DepthTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_DP_SK,Number=1,Type=Float,Description="Jacquard depth for Strelka (based on DP2),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)
 
    def format(self, alt, filter, info, format_dict, count):
        if "DP2" in format_dict.keys():
            format_dict["JQ_DP_SK"] = format_dict["DP2"]
        elif "AU" in format_dict.keys():
            tags = ["AU", "CU", "TU", "GU"]
            denominator = 0
            for tag in tags:
                denominator += int(format_dict[tag].split(",")[1])
            format_dict["JQ_DP_SK"] = str(denominator)
             
        return format_dict

class Strelka_SomaticTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_SOM_SK,Number=1,Type=Integer,Description="Jacquard somatic status for Strelka: 0=non-somatic,1= somatic (based on PASS in FILTER column),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)
 
    def format(self, alt, filter, info, format_dict, count):
        if filter == "PASS":
            format_dict["JQ_HC_SOM_SK"] = self.somatic_status(count)
        else:
            format_dict["JQ_HC_SOM_SK"] = "0"
             
        return format_dict
        
    def somatic_status(self, count):
        if count == 0: #it's NORMAL
            return "0"
        else: #it's TUMOR
            return "1"
        
class LineProcessor():
    def __init__(self, tags):
        self.tags = tags
        
    def add_tags(self, input_line):
        no_newline_line = input_line.rstrip("\n")
        original_vcf_fields = no_newline_line.split("\t")
        new_vcf_fields = original_vcf_fields[:8]
        alt = original_vcf_fields[4]
        filter = original_vcf_fields[6]
        info = original_vcf_fields[7]
        format = original_vcf_fields[8]
        samples = original_vcf_fields[9:]

        count = 0
        for sample in samples:
            format_dict = OrderedDict(zip(format.split(":"), sample.split(":")))
            for tag in self.tags:
                format_dict = tag.format(alt, filter, info, format_dict, count)
                
            if count < 1: ##only add format column once
                new_vcf_fields.append(":".join(format_dict.keys()))
            new_vcf_fields.append(":".join(format_dict.values()))
            count += 1

        return "\t".join(new_vcf_fields) + "\n"

class FileProcessor():
    def __init__(self, tags=[], execution_context_metadataheaders = []):
        self._tags = tags
        self._metaheader = self._metaheader_handler(execution_context_metadataheaders)
        
        for tag in self._tags:
            self._metaheader += tag.metaheader
        self._lineProcessor = LineProcessor(self._tags)
            
    def _metaheader_handler(self, metaheaders):
        new_headers = ["{}\n".format(header) for header in metaheaders]
        return ''.join(new_headers)

    def process(self, reader, writer, caller):
        for line in reader:
            if line.startswith("##"):
                writer.write(line)
            elif line.startswith("#"):
                if caller == "VarScan":
                    if "NORMAL" in line and "TUMOR" in line:
                        writer.write(self._metaheader)
                        writer.write("##jacquard.tag.caller={0}\n".format(caller))
                        writer.write(line)
                    else:
                        print "Unexpected VarScan VCF structure - missing NORMAL\t and TUMOR\n headers."
                        exit(1)
                else:
                    writer.write(self._metaheader)
                    writer.write("##jacquard.tag.caller={0}\n".format(caller))
                    writer.write(line)
            else:
                edited_line = self._lineProcessor.add_tags(line)
                writer.write(edited_line)


def add_subparser(subparser):
    parser_tag = subparser.add_parser("tag", help="Accepts a directory of VCf results and creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser_tag.add_argument("input_dir", help="Path to directory containing VCFs. Other file types ignored")
    parser_tag.add_argument("output_dir", help="Path to Jacquard-tagged VCFs. Will create if doesn't exist and will overwrite files in output directory as necessary")

def determine_file_types(input_dir, in_files, callers):
    file_types = defaultdict(list)   
    inferred_callers = []                                    
    for file in in_files:
        for caller in callers:
            in_file = open(os.path.join(input_dir, file), "r")
            caller_name, valid = caller.validate_input_file(in_file)
            if valid == 1:
                file_types[caller_name].append(file)
                if caller_name == "Unknown":
                    print "ERROR: {0}: ##jacquard.tag.handler={1}".format(os.path.basename(file), caller_name)
                else:
                    print "{0}: ##jacquard.tag.handler={1}".format(os.path.basename(file), caller_name)
                    inferred_caller = "##jacquard.tag.caller={0}".format(caller_name)
                    inferred_callers.append(inferred_caller)
                    print "{0}: {1}".format(os.path.basename(file), inferred_caller)
                break
    return file_types, inferred_callers
    
def print_file_types(file_types):
    for key, vals in file_types.items():
        print "Recognized [{0}] {1} file(s)".format(len(vals), key)
    if "Unknown" in file_types.keys():
        print "Error: unable to determine which caller was used on [{0}]. Check input files and try again.".format(file_types["Unknown"])
        exit(1)
    
def tag_files(input_dir, output_dir, callers, execution_context=[]):
    in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
    if len(in_files) < 1:
        print "Error: Specified input directory [{0}] contains no VCF files. Check parameters and try again."
        exit(1)

    print "\n".join(execution_context)
    print "Processing [{0}] VCF file(s) from [{1}]".format(len(in_files), input_dir)
    
    file_types, inferred_callers = determine_file_types(input_dir, in_files, callers)
    print_file_types(file_types)

    processors = {"VarScan" : FileProcessor(tags=[Varscan_AlleleFreqTag(), Varscan_DepthTag(), Varscan_SomaticTag()], execution_context_metadataheaders=execution_context), "MuTect": FileProcessor(tags=[Mutect_AlleleFreqTag(), Mutect_DepthTag(), Mutect_SomaticTag()], execution_context_metadataheaders=execution_context), "Strelka": FileProcessor(tags=[Strelka_AlleleFreqTag(), Strelka_DepthTag(), Strelka_SomaticTag()], execution_context_metadataheaders=execution_context)}
    
    total_number_of_files = len(in_files)
    count = 1
    for file in in_files:
        print "Reading [{0}] ({1}/{2})".format(os.path.basename(file), count, total_number_of_files)
        fname, extension = os.path.splitext(os.path.basename(file))
        new_file = fname + ".jacquardTags" + extension
        
        in_file = open(os.path.join(input_dir, file), "r")
        out_file = open(os.path.join(output_dir, new_file), "w")
        
        for key, vals in file_types.items():
            if file in vals:
                processors[key].process(in_file, out_file, key)
        
        in_file.close()
        out_file.close()
        
        count += 1
    
    out_files = sorted(glob.glob(os.path.join(output_dir,"*.vcf")))
    
    print "Wrote [{0}] VCF file(s) to [{1}]".format(len(out_files), output_dir)

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    jacquard_utils.validate_directories(input_dir, output_dir)
    
    callers = [Mutect(), Varscan(), Strelka(), Unknown()]
    tag_files(input_dir, output_dir, callers, execution_context)
    
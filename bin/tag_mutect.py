#!/usr/bin/python2.7
import os
from os import listdir
import shutil

class AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_AF_MT,Number=1,Type=Float, Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA).">\n'

    def format(self, format_param_string, format_value_string):
        format_param_array = format_param_string.split(":")
        format_value_array = format_value_string.split(":")
        format_dict = dict(zip(format_param_array, format_value_array))
        
        if "FA" in format_dict.keys():
            format_value_string += ":" + self.roundTwoDigits(format_dict["FA"])
            format_param_string += ":JQ_AF_MT"
            
        return format_param_string, format_value_string

    def roundTwoDigits(self, value): 
        if len(value.split(".")[1]) <= 2:
            return value
        else:
            return str(round(100 * float(value))/100) 
        
class DepthTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_DP_MT,Number=1,Type=Float, Description="Jacquard depth for MuTect (based on DP).">\n'

    def format(self, format_param_string, format_value_string):
        format_param_array = format_param_string.split(":")
        format_value_array = format_value_string.split(":")
        format_dict = dict(zip(format_param_array, format_value_array))
        
        if "DP" in format_dict.keys():
            format_value_string += ":" + format_dict["DP"]
            format_param_string += ":JQ_DP_MT"
            
        return format_param_string, format_value_string
    
class SomaticTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_SOM_MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag).">\n'

    def format(self, format_param_string, format_value_string):
        format_param_array = format_param_string.split(":")
        format_value_array = format_value_string.split(":")
        format_dict = dict(zip(format_param_array, format_value_array))
        
        if "SS" in format_dict.keys():
            format_value_string += ":" + self.somatic_status(format_dict["SS"])
            format_param_string += ":JQ_SOM_MT"
            
        return format_param_string, format_value_string

    def somatic_status(self, ss_value):
        if ss_value == "2":
            return "1"
        else:
            return "0"

class LineProcessor():
    def __init__(self, tags):
        self.tags = tags

    def add_tags(self, input_line):
        line   = input_line.split("\t")[:8]
        format = input_line.split("\t")[8]
        samples = input_line.split("\t")[9:]         

        for tag in self.tags:
            values = []
            for sample in samples:
                param, value = tag.format(format, sample)
                values.append(value)
                
        line.append(param)
        line.extend(values)

        return "\t".join(line)

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
                
#         writer.close()
        
def tag_mutect_files(input_dir, output_dir):
    af_tag = AlleleFreqTag()
    processor = FileProcessor(tags=[af_tag])
    
    for file in sorted(listdir(input_dir)):
        fname, extension = os.path.splitext(file)
        new_file = fname + "_jacquard" + extension
        
        in_file = open(os.path.join(input_dir, file), "r")
        
        for line in in_file:
            if line.startswith("##MuTect"):
                in_file.close()
                in_file = open(os.path.join(input_dir, file), "r")
                out_file = open(os.path.join(output_dir, new_file), "w")
                processor.process(in_file, out_file)
                out_file.close()
                break
        in_file.close()


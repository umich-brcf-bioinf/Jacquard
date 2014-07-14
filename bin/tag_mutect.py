#!/usr/bin/python2.7

class AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_AF_MT,Number=1,Type=Float, Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA).">'

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
        self.metaheader = '##FORMAT=<ID=JQ_DP_MT,Number=1,Type=Float, Description="Jacquard depth for MuTect (based on DP).">'

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
        self.metaheader = '##FORMAT=<ID=JQ_SOM_MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag).">'

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
    def __init__(self, reader, execution_context_metadataheader = []):
        self.reader = reader
        self.metaheader = self._metaheader_handler(execution_context_metadataheader)
            
    def _metaheader_handler(self, metaheaders):
        new_headers = ["##{}\n".format(header) for header in metaheaders]
        return ''.join(new_headers)

    def process(self, writer):
        writer.write(self.metaheader)
        
        for line in self.reader:
            if line.startswith("##"):
                writer.write(line)
        writer.close()
        

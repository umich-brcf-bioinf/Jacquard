#!/usr/bin/python2.7

class AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_AF_MT,Number=1,Type=Float, Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA).">'
        
    def format(self, format_param_string, format_value_string):
        format_param_array = format_param_string.split(":")
        format_value_array = format_value_string.split(":")
        
        format_dict = dict(zip(format_param_array, format_value_array))
        
        final_param_string = format_param_string + ":JQ_AF_MT"
        final_value_string = format_value_string + ":" + format_dict["FA"]
        
        return final_param_string, final_value_string
                
        
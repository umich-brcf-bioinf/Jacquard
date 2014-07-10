#!/usr/bin/python2.7

class AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID=JQ_AF_MT,Number=1,Type=Float, Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA).">'
        
    def format(self, format_param_string, format_value_string):
        format_param_array = format_param_string.split(":")
        format_value_array = format_value_string.split(":")
        format_dict = dict(zip(format_param_array, format_value_array))

        if "FA" in format_dict.keys():
            rounded_FA = "{0:.2f}".format(float(format_dict["FA"])) if len(format_dict["FA"].split(".")[1]) > 2 else format_dict["FA"]
            final_value_string = format_value_string + ":" + rounded_FA
        else:
            final_value_string = format_value_string
            
        final_param_string = format_param_string + ":JQ_AF_MT" if "FA" in format_dict.keys() else format_param_string
       
       
        return final_param_string, final_value_string
                
        
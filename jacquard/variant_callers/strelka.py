#!/usr/bin/python2.7
from collections import defaultdict
import glob
import os
from os import listdir
import re
import jacquard_utils

class Strelka():
    def __init__(self):
        self.name = "Strelka"
        self.meta_header = "##jacquard.normalize_strelka.sources={0},{1}\n"
        self.file_name_search = "snvs|indels"
    
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
    
    def final_steps(self, hc_candidates, merge_candidates, output_dir):
        print "Wrote [{0}] VCF files to [{1}]".format(len(merge_candidates.keys()), output_dir)
        
        return merge_candidates
    
    def handle_hc_files(self, in_file, out_dir, hc_candidates):
        return hc_candidates
        
    def validate_file_set(self, all_keys):
        pass

class AlleleFreqTag():
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

class DepthTag():
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

class SomaticTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}SK,Number=1,Type=Integer,Description="Jacquard somatic status for Strelka: 0=non-somatic,1= somatic (based on PASS in FILTER column),Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_filter_tag, jacquard_utils.__version__)
 
    def format(self, alt, filter, info, format_dict, count):
        strelka_tag = jacquard_utils.jq_filter_tag + "SK"
        if filter == "PASS":
            format_dict[strelka_tag] = self.somatic_status(count)
        else:
            format_dict[strelka_tag] = "0"
             
        return format_dict
        
    def somatic_status(self, count):
        if count == 0: #it's NORMAL
            return "0"
        else: #it's TUMOR
            return "1"
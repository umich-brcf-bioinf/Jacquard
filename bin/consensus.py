#!/usr/bin/python2.7
import argparse
from collections import OrderedDict
import os
from os import listdir
from os.path import isfile, join
import jacquard_utils

def add_consensus(reader, writer, output_file):
    meta_headers = []
    header = ""
    lines = []
    for line in reader:
        if line.startswith("##"):
            meta_headers.append(line)
        elif line.startswith("#"):
            header = line
        else:
            new_line = process_line(line)
            lines.append(new_line)
    
    consensus_meta_headers = ['##FORMAT=<ID=JQ_SOM_SUM,Number=1,Type=Integer,Description="Jacquard consensus somatic call = sum(JQ_SOM_MT, JQ_SOM_SK, JQ_SOM_VS)">', 
                              '##FORMAT=<ID=JQ_AF_AVERAGE,Number=A,Type=Integer,Description="Jacquard consensus somatic call = average(JQ_AF_MT, JQ_AF_SK, JQ_AF_VS)">',
                              '##FORMAT=<ID=JQ_DP_AVERAGE,Number=1,Type=Integer,Description="Jacquard consensus depth = average(JQ_DP_MT, JQ_DP_SK, JQ_DP_VS)">']
    print "\n".join(consensus_meta_headers)
    meta_headers.extend(consensus_meta_headers)
    
    meta_headers.append(header)
    jacquard_utils.write_output(writer, meta_headers, lines)
    print "Wrote consensus-somatic-tagged VCF to [{0}]".format(output_file)
    
def roundTwoDigits(value): 
    new_values = []
    for val in value:
        if len(val.split(".")[1]) <= 2:
            new_values.append(val)
        else:
            new_values.append(str(round(100 * float(val))/100))
    return ",".join(new_values) 
        
def create_consensus_dict(key, val, input_dict, consensus_dict, type):
    split_som = val.split(",")
    for index, af in enumerate(split_som):
        if index in consensus_dict.keys():
            consensus_dict[index] += int(af) if type == "int" else float(af)
        else:
            consensus_dict[index] = int(af) if type == "int" else float(af)
    
    return consensus_dict

def get_consensus_som(field_dict):
    field_list = []
    for key in field_dict.keys():
        field_list.append(str(field_dict[key]))
    consensus = ",".join(field_list) if field_list != [] else 0
    
    return consensus

def get_consensus(consensus_tags, consensus_dict):
    if len(consensus_tags) != 0:
        average = []
        for key in consensus_dict.keys():
            avg = float(consensus_dict[key])/len(consensus_tags)
            average.append(str(avg))
        consensus = ",".join(average)
    else:
        consensus = 0
    
    return consensus

def calculate_consensus(combined_dict):
    consensus_af_tags = []
    consensus_dp_tags = []
    af = {}
    somatic = {}
    depth = {}
    for key in combined_dict.keys():
        if key.startswith("JQ_SOM"):
            if combined_dict[key] != ".":
                somatic = create_consensus_dict(key, combined_dict[key], combined_dict, somatic, "int")
        elif key.startswith("JQ_AF"):
            if combined_dict[key] != ".":
                new_af = roundTwoDigits(combined_dict[key].split(","))
                af = create_consensus_dict(key, new_af, combined_dict, af, "float")
                consensus_af_tags.append(key)
        elif key.startswith("JQ_DP"):
            if combined_dict[key] != ".":
                depth = create_consensus_dict(key, combined_dict[key], combined_dict, depth, "float")
                consensus_dp_tags.append(key)
                
    consensus_af = get_consensus(consensus_af_tags, af)
    consensus_som = get_consensus_som(somatic)
    consensus_dp = get_consensus(consensus_dp_tags, depth)

    combined_dict["JQ_AF_AVERAGE"] = str(consensus_af)
    combined_dict["JQ_SOM_SUM"] = str(consensus_som)
    combined_dict["JQ_DP_AVERAGE"] = str(consensus_dp)
    
    return combined_dict

def process_line(line):
    split_line = line.split("\t")
    format = split_line[8]
    samples = split_line[9:]
    new_samples = []
    for sample in samples:
        combined_dict = jacquard_utils.combine_format_values(format, sample)
        combined_dict = calculate_consensus(combined_dict)
        new_samples.append(":".join(combined_dict.values()))
    new_format = [":".join(combined_dict.keys())]
    new_line = "\t".join(split_line[:8] + new_format + new_samples) + "\n"

    return new_line

def add_subparser(subparser):
    parser_tag = subparser.add_parser("consensus", help="Accepts a Jacquard-merged VCf file and creates a new file, adding consensus fields.")
    parser_tag.add_argument("input_file", help="Path to Jacquard-merged VCF (or any VCF with Jacquard tags (e.g. JQ_SOM_MT)")
    parser_tag.add_argument("output_file", help="Path to output VCf")

def execute(args, execution_context): 
    input_file = os.path.abspath(args.input_file)
    output_file = os.path.abspath(args.output_file)
    print "\n".join(execution_context)
    
    input_file_reader = open(input_file, "r")
    output_file_reader = open(output_file, "w")
    add_consensus(input_file_reader, output_file_reader, output_file)
    
    input_file_reader.close()
    output_file_reader.close()
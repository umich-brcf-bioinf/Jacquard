#!/usr/bin/python2.7
import argparse
from collections import defaultdict, OrderedDict
import numpy
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import re
import jacquard_utils

def calculate_zscore(af_range, dp_range, combined_dict):
    af_mean = sum(af_range)/len(af_range)
    af_std = numpy.std(af_range)
    af_range = float(combined_dict["JQ_AF_RANGE"])
    af_zscore = (af_range - af_mean)/af_std
    rounded_af_zscore = roundTwoDigits([str(af_zscore)])
    
    dp_mean = sum(dp_range)/len(dp_range)
    dp_std = numpy.std(dp_range)
    dp_range = float(combined_dict["JQ_DP_RANGE"])
    dp_zscore = (dp_range - dp_mean)/dp_std
    rounded_dp_zscore = roundTwoDigits([str(dp_zscore)])
    
    combined_dict["JQ_AF_RANGE_ZSCORE"] = str(rounded_af_zscore)
    combined_dict["JQ_DP_RANGE_ZSCORE"] = str(rounded_dp_zscore)
    
    return combined_dict

def iterate_file(reader, writer, output_file, af_range, dp_range, type):
    meta_headers = []
    header = ""
    lines = []
    for line in reader:
        if line.startswith("##"):
            meta_headers.append(line)
        elif line.startswith("#"):
            header = line
        else:
            new_line = process_line(line, af_range, dp_range, type)
            lines.append(new_line)
    
    return meta_headers, header, lines
            
def add_zscore(meta_headers, header, lines, writer, output_file, af_range, dp_range):
    rounded_mean_af = roundTwoDigits([str(sum(af_range)/len(af_range))])
    rounded_mean_dp = roundTwoDigits([str(sum(dp_range)/len(dp_range))])
    rounded_std_af = roundTwoDigits([str(numpy.std(af_range))])
    rounded_std_dp = roundTwoDigits([str(numpy.std(dp_range))])
    
    consensus_meta_headers = ['##FORMAT=<ID= JQ_AF_RANGE_ZSCORE,Number=A,Type=Integer,Description="Jacquard measure of consistency of allele frequencies among callers = (sample AF range - population mean AF range)/standard dev(population AF range)">\n',
                              '##jacquard.consensus.JQ_AF_RANGE_ZSCORE.mean_AF_range={0}\n'.format(rounded_mean_af),
                              '##jacquard.consensus.JQ_AF_RANGE_ZSCORE.standard_deviation={0}\n'.format(rounded_std_af),
                              '##FORMAT=<ID= JQ_DP_RANGE_ZSCORE,Number=A,Type=Integer,Description="Jacquard measure of consistency of depth among callers = (sample DP range - population mean DP range)/standard dev(population DP range)">\n',
                              '##jacquard.consensus.JQ_DP_RANGE_ZSCORE.mean_DP_range={0}\n'.format(rounded_mean_dp),
                              '##jacquard.consensus.JQ_DP_RANGE_ZSCORE.standard deviation_DP_range={0}\n'.format(rounded_std_dp)]
    meta_headers.extend(consensus_meta_headers)
    print "".join(consensus_meta_headers)
    
    meta_headers.append(header)
    jacquard_utils.write_output(writer, meta_headers, lines)
    print "Wrote consensus-somatic-tagged VCF to [{0}]".format(output_file)
    
def add_consensus(meta_headers, header, lines, writer, output_file):
    consensus_meta_headers = ['##FORMAT=<ID=JQ_SOM_SUM,Number=1,Type=Integer,Description="Jacquard consensus somatic call = sum(JQ_SOM_MT, JQ_SOM_SK, JQ_SOM_VS)">\n', 
                              '##FORMAT=<ID=JQ_AF_AVERAGE,Number=A,Type=Integer,Description="Jacquard consensus somatic call = average(JQ_AF_MT, JQ_AF_SK, JQ_AF_VS)">\n',
                              '##FORMAT=<ID=JQ_DP_AVERAGE,Number=1,Type=Integer,Description="Jacquard consensus depth = average(JQ_DP_MT, JQ_DP_SK, JQ_DP_VS)">\n']
    print "".join(consensus_meta_headers)
    meta_headers.extend(consensus_meta_headers)
            
    meta_headers.append(header)
    jacquard_utils.write_output(writer, meta_headers, lines)
    
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

def get_range(consensus_tags, combined_dict, range):
    values = []
    if len(consensus_tags) > 1:
        for tag in consensus_tags:
            values.append(combined_dict[tag])
    if values != []:
        this_range = float(max(values)) - float(min(values))
        range.append(this_range)
    else:
        this_range = 0
        range.append(this_range)

    return range, this_range

def calculate_consensus(combined_dict, af_range, dp_range):
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
    
    af_range, this_af_range = get_range(consensus_af_tags, combined_dict, af_range)
    dp_range, this_dp_range = get_range(consensus_dp_tags, combined_dict, dp_range)
    
    combined_dict["JQ_AF_RANGE"] = str(this_af_range)
    combined_dict["JQ_DP_RANGE"] = str(this_dp_range)
    
    return combined_dict, af_range, dp_range

def combine_format_values(format, sample, sample_column_name):
    new_format = [x + "_" + sample_column_name for x in format.split(":")]
    return OrderedDict(zip(new_format, sample.split(":")))

def process_line(line, af_range, dp_range, type):
    split_line = line.split("\t")
    format = split_line[8]
    samples = split_line[9:]
    new_samples = []

    for sample in samples:
        combined_dict = jacquard_utils.combine_format_values(format, sample)
        
        if type == "zscore":
            combined_dict = calculate_zscore(af_range, dp_range, combined_dict)
        elif type == "consensus":
            combined_dict, af_range, dp_range = calculate_consensus(combined_dict, af_range, dp_range)
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
    
    fname, extension = os.path.splitext(os.path.basename(input_file))
    if not os.path.isfile(input_file) or extension != ".vcf":
        print "ERROR: Input file [{0}] must be a VCF file.".format(input_file)
        exit(1)
        
    fname, extension = os.path.splitext(os.path.basename(output_file))
    if extension != ".vcf":
        print "ERROR: Output file [{0}] must be a VCF file.".format(output_file)
        exit(1)
        
    print "\n".join(execution_context)
    
    af_range = []
    dp_range = []
    
    tmp_file = output_file + ".tmp"
    input_file_reader = open(input_file, "r")
    tmp_file_writer = open(tmp_file, "w")

    print "Adding consensus values to temporary file [{0}]".format(tmp_file)
    meta_headers, header, lines = iterate_file(input_file_reader, tmp_file_writer, tmp_file, af_range, dp_range, "consensus")
    add_consensus(meta_headers, header, lines, tmp_file_writer, tmp_file)
    
    input_file_reader.close()
    tmp_file_writer.close()
    
    tmp_file_reader = open(tmp_file, "r")
    output_file_writer = open(output_file, "w")
    
    print "Adding z-scores to [{0}]".format(output_file)
    meta_headers, header, lines = iterate_file(tmp_file_reader, output_file_writer, output_file, af_range, dp_range, "zscore")
    add_zscore(meta_headers, header, lines, output_file_writer, output_file, af_range, dp_range)
    
    tmp_file_reader.close()
    output_file_writer.close()
    
    os.remove(tmp_file)
    print "Removed temporary file [{0}]".format(tmp_file)
    
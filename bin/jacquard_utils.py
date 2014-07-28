#!/usr/bin/python2.7
__version__ = 0.1

from collections import OrderedDict
from operator import itemgetter, attrgetter
import os
from os import listdir


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
            
def write_output(writer, headers, variants):
    for line in headers:
        writer.write(line)
    for line in variants:
        writer.write(line)
        
def sort_headers(headers):
    meta_headers = []
    field_header = ""
    for header in headers:
        if header.startswith("##"):
            header = header.replace("\t", "")
            meta_headers.append(header)
        else:
            field_header = header

    meta_headers.append(field_header)

    return meta_headers

def sort_data(all_variants):
    new_variants = []
    for variant in all_variants:
        split_variant = variant.split("\t")
        new_line = change_pos_to_int(split_variant)
        new_variants.append(new_line)
        
    sorted_variants = sorted(new_variants, key=itemgetter(0,1,3,4)) #sort by CHROM, POS, REF, ALT
    
    variants = []
    for variant in sorted_variants:
        new_field = [str(field) for field in variant]
        if "chr" not in new_field[0]:
            new_field[0] = "chr" + new_field[0]
        variants.append("\t".join(new_field))
        
    return variants

def change_pos_to_int(split_line):
    new_line = []
    split_line[0] = split_line[0].strip("chr")
    for field in split_line:
        try:
            new_field = int(field)
            new_line.append(new_field)
        except:
            new_line.append(field)
    return new_line

#!/usr/bin/python2.7
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

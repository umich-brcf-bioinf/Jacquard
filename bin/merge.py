#!/usr/bin/python2.7

import datetime
from operator import itemgetter, attrgetter
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
from pandas import *
import re
import sys 
import time

import jacquard_utils

class PivotError(Exception):
    """Base class for exceptions in this module."""
    pass 


class VariantPivoter():
    MISSING_REQUIRED_COLUMNS_ERROR="The columns of specified dataframe do not " +\
        "contain required columns {0}; review input data."

    def __init__(self, rows, combined_df=pd.DataFrame()):
        self._rows = rows
        self._combined_df = combined_df
 
    def add_file(self, sample_file, header_index, file_names):
        initial_df  = create_initial_df(sample_file, header_index)

        if "#CHROM" in initial_df.columns:
            initial_df.rename(columns={"#CHROM": "CHROM"}, inplace=True)

        unpivoted_df = self.is_compatible(initial_df)
        self._combined_df = merge_samples(unpivoted_df, self._combined_df, self._rows, file_names)

    def validate_sample_data(self):
        grouped = self._combined_df.groupby(self._rows)
        group = grouped.groups

        for key, val in group.items():
            if len(val) != 1:
                for column in self._combined_df:
                    self.find_non_unique_rows(grouped, column, key, val)
        for column in self._combined_df:
            self.find_non_unique_cells(column)
        

        self._combined_df.reset_index(inplace=True)   
        del self._combined_df["index"]
        
        return self._combined_df
        
    def find_non_unique_rows(self, grouped, column, key, val):
        col_data = []
        for index in val:
            data = grouped.get_group(key).ix[index, column]
            if data not in col_data:
                col_data.append(data)
        if len(col_data) != 1:
            self._combined_df.ix[val[-1], column] = "^" 
            self._combined_df = self._combined_df.drop(self._combined_df.index[val[:-1]])
        
    def find_non_unique_cells(self, column):
        count = 0
        for data in self._combined_df[column]:
            if type(data) == np.ndarray:
                self._combined_df.ix[count, column] = "^"
            count += 1
                
    def sort_rows(self, df):
        df.reset_index(inplace=True)
        sorted_df = df.sort(self._rows)

        return sorted_df
                
    def is_compatible(self, initial_df):   
        if "#CHROM" in initial_df.columns:
            initial_df.rename(columns={"#CHROM": "CHROM"}, inplace=True)
            
        initial_df = project_prepivot(initial_df)
        
        self._check_required_columns_present(initial_df)
        self._check_pivot_is_unique(initial_df)  
       
        return initial_df

    def _check_required_columns_present(self, dataframe):
        required_columns = set(self._rows)

        if not required_columns.issubset(dataframe.columns.values):
            raise PivotError("Missing required columns; contact sysadmin.")
         
    def _check_pivot_is_unique(self, dataframe):   
        group = self._rows
        grouped_df = dataframe.groupby(group)

        if len(grouped_df.groups) != len(dataframe):
            raise PivotError("Duplicate keys would result in an invalid pivot; contact sysadmin.")

def validate_parameters(input_keys, first_line, header_names):
    invalid_fields = []
    fields = header_names.split("\t")

    for key in input_keys:
        if key not in fields:
            invalid_fields.append(key)
    
    message = "Invalid input parameter(s) "
    raise_err = 0
    if invalid_fields != []:
        message += str(invalid_fields)
        raise_err = 1
    
    return raise_err, message

def project_prepivot(df):
    remove_columns = set(["ID", "QUAL", "FILTER", "INFO"])

    for col in df.columns:
        if col in remove_columns:
            del df[col]

    return df

def build_pivoter(sample_file, input_keys, header_index):
    initial_df  = create_initial_df(sample_file, header_index)
    
    pivoter = VariantPivoter(input_keys)
    pivoter.is_compatible(initial_df)
    
    return pivoter

def create_initial_df(sample_file, header_index):
    initial_df = pd.read_csv(sample_file, sep="\t", header=header_index, dtype='str')

    return initial_df
  
def merge_samples(reduced_df, combined_df, rows, file_names):
    if combined_df.empty:
        combined_df = reduced_df
    else:
        combined_df = merge(combined_df, reduced_df, suffixes=file_names, how="outer", on=rows)

    return combined_df    
    
def rearrange_columns(output_df):
    index = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"] 
    format = []
    samples = []
    for i in list(output_df.columns.values):
        if i in index:
            continue
        elif re.search("FORMAT", i):
            format.append(i)
        else:
            samples.append(i)
            
    headers = index + format + samples
    ##change order of columns:
    output_df = output_df.ix[:,headers]
    return output_df
    
def process_files(sample_file_readers, input_dir, output_path, input_keys, headers, header_names, first_line, pivot_builder=build_pivoter):
    first_file_reader = sample_file_readers[0]
    first_file      = first_file_reader
    
    print "{0} - validating command line parameters".format(datetime.fromtimestamp(time.time()).strftime('%Y/%m/%d %H:%M:%S'))
    raise_err, message = validate_parameters(input_keys, first_line, header_names)  
    if raise_err == 1:
        raise PivotError(message)
        
    print "{0} - determining file type".format(datetime.fromtimestamp(time.time()).strftime('%Y/%m/%d %H:%M:%S'))
    pivoter  = pivot_builder(first_file, input_keys, headers[0])
    
    file_names = []
    count = 0
    for sample_file in sample_file_readers:
        print "{0} - reading ({1}/{2}): {3}".format(datetime.fromtimestamp(time.time()).strftime('%Y/%m/%d %H:%M:%S'), count + 1, len(sample_file_readers), sample_file)
        
        file_names.append("_" + os.path.basename(sample_file))
        pivoter.add_file(sample_file, headers[count], file_names)
        
        count += 1
   
    print "{0} - validating sample data".format(datetime.fromtimestamp(time.time()).strftime('%Y/%m/%d %H:%M:%S'))
    pivoter.validate_sample_data()

    print "{0} - rearranging columns".format(datetime.fromtimestamp(time.time()).strftime('%Y/%m/%d %H:%M:%S'))
    complete_output = rearrange_columns(pivoter._combined_df)
#     
    print "{0} - sorting rows".format(datetime.fromtimestamp(time.time()).strftime('%Y/%m/%d %H:%M:%S'))
    sorted_df = pivoter.sort_rows(complete_output)
# 
    if "index" in sorted_df:
        del sorted_df["index"]
        
    sorted_df = sorted_df.fillna(".")

    print "{0} - writing to vcf file: {1}".format(datetime.fromtimestamp(time.time()).strftime('%Y/%m/%d %H:%M:%S'), output_path)
    sorted_df.to_csv(output_path, index=False, sep="\t")  

#    
    print "{0} - done".format(datetime.fromtimestamp(time.time()).strftime('%Y/%m/%d %H:%M:%S'))
    
def determine_input_keys(input_dir):
    for file in listdir(input_dir):
        if isfile(join(input_dir, file)):
            fname, extension = os.path.splitext(file)
            if extension == ".vcf":
                return ["CHROM", "POS", "REF", "ALT"]
                
            else:
                raise PivotError("Cannot determine columns to be used as keys for the pivoting from file [{0}]. Please specify parameter [-k] or [--keys]".format(os.path.abspath(file)))
    
def get_headers_and_readers(input_dir):
    sample_file_readers = []
    headers = []
    header_names = []
    first_line = []
    
    for item in sorted(listdir(input_dir)):
        if isfile(join(input_dir, item)):
            f = open(input_dir + "/" + item, 'r')
            count = -1
            for line in f:
                count += 1
                if line.startswith("##"):
                    continue
                elif line.startswith("#"):
                    headers.append(count)
                    header_names.append(line)
                else:
                    first_line.append(line)
                    break

            f.close()
            sample_file_readers.append(input_dir + "/" + item)

    header_names = header_names[0]
    
    header_names = re.sub(r'#CHROM', 'CHROM', header_names)
    
    return sample_file_readers, headers, header_names, first_line

def add_subparser(subparser):
    parser_pivot = subparser.add_parser("merge", help="Accepts a directory of VCFs and returns a single merged VCF file.")
    parser_pivot.add_argument("input_dir", help="Path to directory containing VCFs. Other file types ignored")
    parser_pivot.add_argument("output_file", help="Path to output variant-level VCF file")
    parser_pivot.add_argument("-k", "--keys",
        help="Columns to be used as keys for the pivoting. Default keys for VCF are CHROM,POS,ID,REF,ALT,QUAL,FILTER.")
        
def execute(args, execution_context):
    input_dir = os.path.abspath(args.input_dir)
    output_path = os.path.abspath(args.output_file)
    input_keys = args.keys.split(",") if args.keys else determine_input_keys(input_dir)
    
    output_dir, outfile_name = os.path.split(output_path)

    jacquard_utils.validate_directories(input_dir, output_dir)
        
    fname, extension = os.path.splitext(outfile_name)
    if extension != ".vcf": 
        print "Error. Specified output {0} must have a .vcf extension".format(output_path)
        exit(1)
    
    sample_file_readers, headers, header_names, first_line = get_headers_and_readers(input_dir)
    process_files(sample_file_readers, input_dir, output_path, input_keys, headers, header_names, first_line)
    
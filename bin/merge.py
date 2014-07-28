#!/usr/bin/python2.7
from collections import defaultdict, OrderedDict
import datetime
import glob
import math
from operator import itemgetter, attrgetter
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
from pandas import *
import re
from sets import Set
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
 
    def add_file(self, sample_file, header_index, file_name):
        initial_df  = create_initial_df(sample_file, header_index)

        if "#CHROM" in initial_df.columns:
            initial_df.rename(columns={"#CHROM": "CHROM"}, inplace=True)

        unpivoted_df = self.is_compatible(initial_df)
        fname_df, sample_columns = append_fname_to_samples(initial_df, file_name, self._rows)
        self._combined_df = merge_samples(fname_df, self._combined_df, self._rows)
        
        return sample_columns

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
  
def append_fname_to_samples(df, file_name, rows):
    indexed_df = df.set_index(rows)
    basename = os.path.basename(file_name)
    new_df = indexed_df.rename(columns=lambda x: "|".join([basename, x]))
    sample_columns = new_df.columns.values
    reset_df = new_df.reset_index()
    
    return reset_df, sample_columns
  
def merge_samples(reduced_df, combined_df, rows):
    if combined_df.empty:
        combined_df = reduced_df
    else:
        combined_df = merge(combined_df, reduced_df, how="outer", on=rows)
        
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
    
def create_dict(df, row, columns):
    file_dict = defaultdict(list)
    all_tags = []
    for column in columns:
        if re.search("\|", column) and not re.search("\|FORMAT", column):
            fname = column.split("|")[0]

            format_column = str(df.ix[row, fname + "|FORMAT"])
            sample_column = str(df.ix[row, column])
            format_column += ":sample_name"
            sample_column += ":" + column
            format_sample = "{0}={1}".format(format_column, sample_column)
            tags, format_sample_dict = combine_format_values(format_sample)

            file_dict[fname].append(format_sample_dict)
            for tag in tags:
                try: 
                    float(tag)
                except:
                    if tag not in all_tags:
                        all_tags.append(tag)
    return file_dict, all_tags

def combine_format_values(aggregate_col):
    pairs = aggregate_col.split("=")
    tags = pairs[0].split(":")
    return tags, OrderedDict(zip(pairs[0].split(":"), pairs[1].split(":")))

def add_all_tags(file_dict, sample_keys):
    for sample_list in file_dict.values():
        for sample in sample_list:
            for samp_key in sample_keys:
                if samp_key not in sample.keys():
                    sample[samp_key] = "."
    return file_dict

def sort_format_tags(file_dict):
    sorted_file_dict = defaultdict(list)
    for fname, sample_dicts in file_dict.items():
        new_samp_dicts = []
        for sample_dict in sample_dicts:
            new_sample_dict = OrderedDict(sorted(sample_dict.iteritems()))
            new_samp_dicts.append(new_sample_dict)
        sorted_file_dict[fname] = new_samp_dicts
#     
    return sorted_file_dict

def remove_non_jq_tags(file_dict):
    sample_keys = []
    for sample_list in file_dict.values():
        for sample in sample_list:
            for key, val in sample.items():
                if not re.search("JQ_", key) and key != "sample_name":
                    del sample[key]
                else:
                    if key not in sample_keys:
                        sample_keys.append(key)
                        
    file_dict = add_all_tags(file_dict, sample_keys)
    sort_format_tags(file_dict)
    return file_dict

def cleanup_df(df, file_dict):
    for key in file_dict.keys():
        try:
            del df[key + "|FORMAT"]
        except:
            pass
    for col in df.columns:
        df = df.applymap(lambda x: str(x).replace(":" + col, ""))
        df = df.applymap(lambda x: str(x).replace(col + ":", ""))
        df = df.applymap(lambda x: str(x).replace(col, "."))
    
    df = df.applymap(lambda x: str(x).replace(":sample_name", ""))
    df = df.applymap(lambda x: str(x).replace("sample_name:", ""))
    df = df.applymap(lambda x: str(x).replace("sample_name", "."))
    
    df.replace("nan", ".", inplace=True)

    return df

def combine_format_columns(df):
    for row, col in df.T.iteritems():
        columns = col.index.values
        file_dict, all_tags = create_dict(df, row, columns)
#         file_dict = get_consensus_format_sample(file_dict, all_tags)
        file_dict = remove_non_jq_tags(file_dict)
        
        for key, val in file_dict.items():
            for thing in val:
                df.ix[row, "FORMAT"] = ":".join(thing.keys())
                df.ix[row, thing["sample_name"]] = ":".join(thing.values())
            
    df = cleanup_df(df, file_dict)

    return df
    
def determine_merge_execution_context(all_merge_context, all_merge_column_context, sample_columns, sample_file, count):
    samples = []
    actual_sample_columns = []
    
    samp_count = 0
    for samp_column in sample_columns:
        samp_name = samp_column.split("|")[1]
        if samp_name != "FORMAT":
            samp_count += 1
            samples.append(samp_name)
            actual_sample_columns.append(samp_column)
            merge_column_context = "##jacquard.merge.sample_column{0}={1}".format(samp_count, samp_column)
            all_merge_column_context.append(merge_column_context)
            
    merge_context = "##jacquard.merge.file{0}={1}({2})".format(count, os.path.basename(sample_file), samples)
    all_merge_context.append(merge_context)

    return all_merge_context, all_merge_column_context
    
def print_new_execution_context(execution_context, out_file):
    execution_context[-1] = execution_context[-1] + "\n"

    out_file.write("\n".join(execution_context))
    out_file.close()
    
def process_files(sample_file_readers, input_dir, output_path, input_keys, headers, header_names, first_line, execution_context, pivot_builder=build_pivoter):
    first_file_reader = sample_file_readers[0]
    first_file      = first_file_reader
    
    raise_err, message = validate_parameters(input_keys, first_line, header_names)  
    if raise_err == 1:
        raise PivotError(message)
        
    pivoter  = pivot_builder(first_file, input_keys, headers[0])
    
    print "Processing [{0}] VCF files from [{1}]".format(len(sample_file_readers), input_dir)
    count = 0
    all_merge_context = []
    all_merge_column_context = []
    for sample_file in sample_file_readers:
#         print "{0} - reading ({1}/{2}): {3}".format(datetime.fromtimestamp(time.time()).strftime('%Y/%m/%d %H:%M:%S'), count + 1, len(sample_file_readers), sample_file)
        sample_columns = pivoter.add_file(sample_file, headers[count], sample_file)
        count += 1
        
        all_merge_context, all_merge_column_context = determine_merge_execution_context(all_merge_context, all_merge_column_context, sample_columns, sample_file, count)
    
    new_execution_context = all_merge_context + all_merge_column_context
    print "\n".join(new_execution_context)
    
    execution_context.extend(new_execution_context)
    writer = open(output_path, "w")
    print_new_execution_context(execution_context, writer)
# 
    pivoter.validate_sample_data()
    formatted_df = combine_format_columns(pivoter._combined_df)
    rearranged_df = rearrange_columns(formatted_df)
    sorted_df = pivoter.sort_rows(rearranged_df)

    if "index" in sorted_df:
        del sorted_df["index"]
        
    sorted_df = sorted_df.fillna(".")
    
    with open(output_path, "a") as f:
        sorted_df.to_csv(f, index=False, sep="\t")  

    print "Merged [{0}] VCf files to [{1}]".format(len(sample_file_readers), output_path)
    
def determine_input_keys(input_dir):
    for file in listdir(input_dir):
        if isfile(join(input_dir, file)):
            fname, extension = os.path.splitext(file)
            if extension == ".vcf":
                return ["CHROM", "POS", "REF", "ALT"]
                
            else:
                raise PivotError("Cannot determine columns to be used as keys for the pivoting from file [{0}]. Please specify parameter [-k] or [--keys]".format(os.path.abspath(file)))
    
def get_headers_and_readers(in_files):
    sample_file_readers = []
    headers = []
    header_names = []
    first_line = []
    meta_headers = []
    invalid_files = []
    for in_file in in_files:
        f = open(in_file, 'r')
        count = -1
        invalid = 1
        for line in f:
            count += 1
            if line.startswith("##"):
                if re.search("##FORMAT=<ID=JQ_", line):
                    meta_headers.append(line.strip("\n"))
                if re.search("##jacquard.tag.caller=", line):
                    invalid = 0
            elif line.startswith("#"):
                headers.append(count)
                header_names.append(line)
            else:
                first_line.append(line)
                break
        if invalid == 1:
            invalid_files.append(file)

        f.close()
        sample_file_readers.append(in_file)

    if invalid_files != []:
        print "ERROR: VCF file(s) [{0}] have no Jacquard tags. Run [jacard tag] on these files and try again.".format(invalid_files)
        exit(1)

    header_names = header_names[0]
    header_names = re.sub(r'#CHROM', 'CHROM', header_names)
    
    return sample_file_readers, headers, header_names, first_line, meta_headers

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
        
    in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
    if len(in_files) < 1:
        print "Error: Specified input directory [{0}] contains no VCF files. Check parameters and try again."
        exit(1)
        
    sample_file_readers, headers, header_names, first_line, meta_headers = get_headers_and_readers(in_files)
            
    print "\n".join(execution_context) 
    execution_context.extend(meta_headers + ["##fileformat=VCFv4.2"])
    process_files(sample_file_readers, input_dir, output_path, input_keys, headers, header_names, first_line, execution_context)
    
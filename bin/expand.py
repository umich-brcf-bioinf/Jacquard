#!/usr/bin/python2.7
import argparse
import datetime
import numpy
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
from pandas import *
import random
import re
import sys 
import time
import openpyxl
from openpyxl import load_workbook
from openpyxl.style import Color, Fill, Font

import jacquard_utils

class PivotError(Exception):
    """Base class for exceptions in this module."""
    pass 


class VariantPivoter():
    MISSING_REQUIRED_COLUMNS_ERROR="The columns of specified dataframe do not " +\
        "contain required columns {0}; review input data."

    @staticmethod
    def _build_transform_method(rows, columns, pivot_values):
        def transform(df, fname):
            expanded_df = expand_format(df, pivot_values, rows, fname)
            
            return project_prepivot(expanded_df, pivot_values, rows, columns)
        return transform

    def __init__(self, rows, cols, pivot_values, combined_df=pd.DataFrame(), annot_df=pd.DataFrame()):
        self._rows = rows
        self._cols = cols
        self._pivot_values = pivot_values
        self._transform = VariantPivoter._build_transform_method(self._rows, self._cols, self._pivot_values)
        self._combined_df = combined_df
        self._annot_df = annot_df
        
    def pivot(self):
        pivoted_df = pd.pivot_table(
            self._combined_df, values=self._pivot_values, rows=self._rows, cols=self._cols, 
            aggfunc=lambda x: x)

        try:
            pivoted_df = pivoted_df.applymap(lambda x: int(x))
        except:
            pass

        return pivoted_df
    
    def add_file(self, input_file, header_index):
        initial_df  = create_initial_df(input_file, header_index)

        if "#CHROM" in initial_df.columns:
            initial_df.rename(columns={"#CHROM": "CHROM"}, inplace=True)
        
        self._annot_df = append_to_annot_df(initial_df, self._annot_df)
        
    def validate_annotations(self):
        grouped_df = self._annot_df.groupby(self._rows)
        self._annot_df = grouped_df.last()
        self._annot_df.reset_index(inplace=True)
        
        return self._annot_df
                
    def join_dataframes(self, pivot_df):
        output = merge(self._annot_df, pivot_df, left_on=self._rows, right_index=True, how='outer', sort=False)
        output = output.fillna("")

        return output
    
    def sort_rows(self, df):
        df.reset_index(inplace=True)
        sorted_df = df.sort(self._rows)

        return sorted_df
                
    def is_compatible(self, initial_df, fname=""):   
        unpivoted_df = self._transform(initial_df, fname)

        self._check_required_columns_present(unpivoted_df)
        self._check_pivot_is_unique(unpivoted_df)  
       
        return unpivoted_df

    def _check_required_columns_present(self, dataframe):
        required_columns = set(self._rows + self._cols)
        if not required_columns.issubset(dataframe.columns.values):
            print "columns absent"
            raise PivotError("Missing required columns; contact sysadmin.")
         
    def _check_pivot_is_unique(self, dataframe):   
        group = self._rows + self._cols
        grouped_df = dataframe.groupby(group)

        if len(grouped_df.groups) != len(dataframe):
            print "not unique"
            raise PivotError("Duplicate keys would result in an invalid pivot; contact sysadmin.")

def validate_parameters(input_keys, first_line, header_names, pivot_values):
    invalid_fields = []

    fields = header_names.split("\t")

    for key in input_keys:
        if key not in fields:
            if key == "SnpEff_WARNING/ERROR":
                input_keys.remove(key)
                input_keys.append("WARNING/ERROR")
            else:
                invalid_fields.append(key)
    
    invalid_tags = validate_format_tags(first_line, pivot_values, fields)
    
    message = "Invalid input parameter(s) "
    raise_err = 0
    if invalid_fields != []:
        message += str(invalid_fields)
        raise_err = 1
    if invalid_tags != []:
        message += str(invalid_tags)
        raise_err = 1
    
    return raise_err, message

def validate_format_tags(first_line, pivot_values, fields):
    invalid_tags = []
    all_format_tags ={}
    
    for line in first_line:
        my_line = line.split("\t")
        for index, field in enumerate(fields):
            if field == "FORMAT":
                format = my_line[index]
                format_tags = format.split(":")
                
                for tag in format_tags:
                    all_format_tags[tag] = 0

    for val in pivot_values:
        if val not in all_format_tags:
            invalid_tags.append(val)

    return invalid_tags

def create_initial_df(input_file, header_index):
    initial_df = pd.read_csv(input_file, sep="\t", header=header_index, dtype='str')

    return initial_df

def build_pivoter(input_file, input_keys, pivot_values, header_index):
    initial_df  = create_initial_df(input_file, header_index)

    if "SnpEff_WARNING/ERROR" in initial_df.columns.values or "WARNING/ERROR" in initial_df.columns.values:
        initial_df.rename(columns={"WARNING/ERROR": "SnpEff_WARNING/ERROR"}, inplace=True)
        initial_df = _exclude_errors_and_warnings(initial_df)
    
    pivoter = VariantPivoter(input_keys, ["SAMPLE_NAME"], pivot_values)
    unpivoted_df = pivoter.is_compatible(initial_df)
    pivoter._combined_df = unpivoted_df
    
    return pivoter
  
def _exclude_errors_and_warnings(df):
    try:
        filtered_df = df[df['SnpEff_WARNING/ERROR']=='.']
        return filtered_df
    except KeyError as e:
        raise PivotError("File is missing SnpEff_WARNING/ERROR column; review input files")
    else:
        raise
    
def append_to_annot_df(initial_df, annot_df):  
    format_index = initial_df.columns.get_loc("FORMAT")
    last_column = initial_df.shape[1]
    
    annotation_columns = []
    annotation_columns.extend(initial_df.columns)
    
    ##remove column if it is format column or any column to the right of format
    for i in range(format_index, last_column):
        annotation_columns.remove(initial_df.columns[i])
    
    temp_df = pd.DataFrame()
    for col in initial_df.columns:
        temp_df[col] = initial_df[col]
        if col not in annotation_columns:
            del temp_df[col]
            
    if annot_df.empty:
        annot_df = temp_df
    else:
        annot_df = annot_df.append(temp_df, ignore_index=True)
    return annot_df
  
def melt_samples(df, fname): 
    format_data = df.ix[0, "FORMAT"].split(":")

    first_sample_column_index = df.columns.get_loc("FORMAT") + 1
    first_sample_column = df.ix[:, first_sample_column_index]
    last_column = df.shape[1]
    
    sample_names = []
    for i in range(first_sample_column_index, last_column):
        field_name = df.columns[i]
        data = df.ix[0, field_name].split(":")
        if len(data) == len(format_data):
            sample_names.append(field_name)
            

    column_list = []
    column_list.extend(df.columns)

    for col in sample_names:
        if col in column_list: 
            column_list.remove(col)   
        df.rename(columns={col: fname + "_" + col}, inplace=True)
        
    melted_df = melt(df, id_vars=column_list, var_name='SAMPLE_NAME', value_name='SAMPLE_DATA')
    
    return melted_df
    
def expand_format(df, formats_to_expand, rows, fname):
    df = melt_samples(df, fname)

    df["aggregate_format_sample"] = df["FORMAT"] + "=" + df['SAMPLE_DATA']
    df["aggregate_format_sample"] = df["aggregate_format_sample"].map(combine_format_values)

    s = df["aggregate_format_sample"].apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = "format_sample"

    unpivoted_format_value_df = s.apply(lambda x: pd.Series(x))
    unpivoted_format_value_df.columns = ["FORMAT2", "VALUE2"]
    
    joined_df = df.join(unpivoted_format_value_df)
    del joined_df["aggregate_format_sample"]
    
    if "#CHROM" in joined_df.columns:
        joined_df.rename(columns={"#CHROM": "CHROM"}, inplace=True)
        
    if "SnpEff_WARNING/ERROR" in joined_df.columns and "SnpEff_WARNING/ERROR" not in rows:
        joined_df.rename(columns={"SnpEff_WARNING/ERROR": "WARNING/ERROR"}, inplace=True)

    try:
        pivoted_df = pd.pivot_table(joined_df, rows=rows+["SAMPLE_NAME"], cols="FORMAT2", values="VALUE2", aggfunc=lambda x: x)
    except Exception as e :
        raise PivotError("Cannot pivot data. {0}".format(e))


    pivoted_df.reset_index(inplace=True)

    return pivoted_df

def combine_format_values(aggregate_col):
    pairs = aggregate_col.split("=")
    return zip(pairs[0].split(":"), pairs[1].split(":"))
    
def project_prepivot(df, pivot_values, rows, columns):
    required_columns = set(pivot_values + rows + columns)

    for col in df.columns:
        if col not in required_columns:
            del df[col]

    return df   
    
def rearrange_columns(output_df):
    ##change tuples to strings:
    lst = []
    pivot_columns = []
    for i in list(output_df.columns.values):
        if type(i) is tuple:
            i = "_".join(i)
            pivot_columns.append(i)
        lst.append(i)
    output_df.columns = lst 
    
    ##change order of columns:
    lst = change_order(lst, pivot_columns)    
    output_df = output_df.ix[:,lst]

    return output_df
    
def change_order(lst, pivot_columns):
    meta_lst = []
    annot_lst = []
    for header in lst:    
        if re.search("CHROM|POS|REF|ALT|ANNOTATED_ALLELE|GENE_SYMBOL|Mult_Alt|Mult_Gene|IGV|UCSC|dbSNP", header):
            meta_lst.append(header)
        elif header not in pivot_columns:
            annot_lst.append(header) 

    lst = meta_lst + pivot_columns + annot_lst

    return lst
    
def insert_links(joined_output):
    for row, col in joined_output.T.iteritems():
        joined_output.loc[row, "UCSC"] = '=hyperlink("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr' + joined_output.loc[row, "CHROM"] + ":" + str(int(joined_output.loc[row, "POS"]) - 50)+ "-" + str(int(joined_output.loc[row, "POS"]) + 50) + '", "UCSC")'
        joined_output.loc[row, "IGV"] = '=hyperlink("localhost:60151/goto?locus=chr' + joined_output.loc[row, "CHROM"] + ':' + joined_output.loc[row, "POS"] + '", "IGV")'
        joined_output.loc[row, "dbSNP"] = '=hyperlink("http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=&' + joined_output.loc[row, "ID"] + '", "dbSNP")' if joined_output.loc[row, "ID"] != "." else ""
    del joined_output["ID"]
        
def process_files(input_file, output_path, input_keys, pivot_values, headers, header_names, first_line, pivot_builder=build_pivoter):
    raise_err, message = validate_parameters(input_keys, first_line, header_names, pivot_values)  
    if raise_err == 1:
        raise PivotError(message)
        
    pivoter  = pivot_builder(input_file, input_keys, pivot_values, headers[0])
    pivoter.add_file(input_file, headers[0])
    pivoter.validate_annotations()

    pivoted_df = pivoter.pivot()
    pivoted_df = pivoted_df.fillna("")

    joined_df = pivoter.join_dataframes(pivoted_df)
    insert_links(joined_df)

    rearranged_df = rearrange_columns(joined_df)
    try:
        rearranged_df.ix[:, "CHROM"] = rearranged_df.ix[:, "CHROM"].apply(lambda x: int(x.strip("chr")))
    except:
        rearranged_df.ix[:, "CHROM"] = rearranged_df.ix[:, "CHROM"].apply(lambda x: x.strip("chr"))
    rearranged_df.ix[:, "POS"] = rearranged_df.ix[:, "POS"].apply(lambda x: int(x))
    
    sorted_df = pivoter.sort_rows(rearranged_df)
    sorted_df.ix[:, "CHROM"] = sorted_df.ix[:, "CHROM"].apply(lambda x: "chr" + str(x))
    sorted_df.ix[:, "POS"] = sorted_df.ix[:, "POS"].apply(lambda x: str(x))

    if "index" in sorted_df:
        del sorted_df["index"]
         
    if "WARNING/ERROR" in sorted_df.columns.values:
        sorted_df.rename(columns={"WARNING/ERROR":"SnpEff_WARNING/ERROR"}, inplace=True)
     
    writer = ExcelWriter(output_path)
    sorted_df.to_excel(writer, "Variant_output", index=False, merge_cells=0)  
    writer.save() 

    print "Wrote formatted Excel file to [{0}]".format(output_path)
   
def determine_input_keys(input_file):
    fname, extension = os.path.splitext(input_file)

    if extension == ".txt":
        return ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL", "SnpEff_WARNING/ERROR"]
   
    elif extension == ".vcf":
        return ["CHROM", "POS", "REF", "ALT"]
        
    else:
        raise PivotError("Cannot determine columns to be used as keys for the pivoting from file [{0}]. Please specify parameter [-k] or [--keys]".format(os.path.abspath(input_file)))

def get_headers(input_file):
    meta_headers = []
    headers = []
    header_names = []
    first_line = []
    
    f = open(input_file, 'r')
    count = -1
    for line in f:
        count += 1
        if line.startswith("##") or line.startswith("#Epee"):
            meta_headers.append(line)
        elif line.startswith("#"):
            headers.append(count)
            header_names.append(line)
        else:
            first_line.append(line)
            break

    f.close()

    header_names = header_names[0]
    header_names = re.sub(r'#CHROM', 'CHROM', header_names)
    
    return meta_headers, headers, header_names, first_line

def add_subparser(subparser):
    parser_pivot = subparser.add_parser("expand", help="Pivots annotated VCF file so that given sample specific information is fielded out into separate columns. Returns an Excel file containing concatenation of all input files.")
    parser_pivot.add_argument("input_file", help="Path to annotated VCF. Other file types ignored")
    parser_pivot.add_argument("output_file", help="Path to output variant-level XLSX file")
    parser_pivot.add_argument("-k", "--keys",
        help="Columns to be used as keys for the pivoting. Default keys for VCF are CHROM,POS,REF,ALT. Default keys for Epee TXT are CHROM,POS,REF,ANNOTATED_ALLELE,GENE_SYMBOL")
    parser_pivot.add_argument("-t", "--tags",
        help="Format tags to be fielded out in the pivoting.")
        
def execute(args, execution_context):
    input_file = os.path.abspath(args.input_file)
    output_path = os.path.abspath(args.output_file)
    input_keys = args.keys.split(",") if args.keys else determine_input_keys(input_file)
    pivot_values = args.tags.split(",") if args.tags else ["GT"]
    
    output_dir, outfile_name = os.path.split(output_path)
    try:
        os.mkdir(output_dir)
    except:
        pass
    
    fname, extension = os.path.splitext(outfile_name)
    if extension != ".xlsx": 
        print "Error. Specified output {0} must have a .xlsx extension".format(output_path)
        exit(1)
    
    meta_headers, headers, header_names, first_line = get_headers(input_file)
    print "\n".join(execution_context) 
    execution_context.extend(meta_headers + ["##fileformat=VCFv4.2"])
    process_files(input_file, output_path, input_keys, pivot_values, headers, header_names, first_line)
    
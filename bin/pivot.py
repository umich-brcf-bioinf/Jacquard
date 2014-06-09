#!/usr/bin/python2.7
import math
import numpy
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
from pandas import *
import random
import re
import sys 
import openpyxl
from openpyxl import load_workbook
from openpyxl.style import Color, Fill, Font

class PivotError(Exception):
    """Base class for exceptions in this module."""
    pass 


class VariantPivoter():

    MISSING_REQUIRED_COLUMNS_ERROR="The columns of specified dataframe do not " +\
        "contain required columns {0}; review input data."
    NOOP_TRANSFORM = lambda x : x

    def __init__(self, rows, cols, pivot_values, transform=NOOP_TRANSFORM, combined_df=pd.DataFrame(), annot_df=pd.DataFrame()):
        self._rows = rows
        self._cols = cols
        self._pivot_values = pivot_values
        self._transform = transform
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
 
    def add_file(self, path, reader, header_index):
        initial_df  = create_initial_df(path, reader, header_index)
        
        if "#CHROM" in initial_df.columns:
            initial_df.rename(columns={"#CHROM": "CHROM"}, inplace=True)
        
        self._annot_df = append_to_annot_df(initial_df, self._annot_df)
        
        unpivoted_df = self._is_compatible(initial_df)

        self._combined_df = merge_samples(unpivoted_df, self._combined_df)
        
    def validate_annotations(self):
        grouped_df = self._annot_df.groupby(self._rows)
        self._annot_df = grouped_df.last()
        self._annot_df.reset_index(inplace=True)
        
    def validate_sample_data(self):
        grouped = self._combined_df.groupby(self._rows + ["SAMPLE_NAME"])
        group = grouped.groups
        
        for column in self._combined_df:
            for key, val in group.items():
                self.find_non_unique_columns(grouped, column, key, val)
            self.find_non_unique_rows(column)
        return self._combined_df
        
    def find_non_unique_columns(self, grouped, column, key, val):
        if len(val) != 1:
            col_data = []
            for index in val:
                data = grouped.get_group(key).ix[index, column]
                if data not in col_data:
                    col_data.append(data)
            if len(col_data) != 1:
                for index in val:
                    self._combined_df.ix[index, column] = "^"    
                    
    def find_non_unique_rows(self, column):
        count = 0
        for data in self._combined_df[column]:
            if type(data) == np.ndarray:
                self._combined_df.ix[count, column] = "^"
            count += 1
                
    def join_dataframes(self, pivot_df):
        output = merge(self._annot_df, pivot_df, left_on=self._rows, right_index=True, how='outer', sort=False)
        output = output.fillna("")

        return output
    
    def sort_rows(self, df):
        df.reset_index(inplace=True)
        sorted_df = df.sort(self._rows)

        return sorted_df
    
    def insert_mult_alt(self, df):
        df.insert(0, "Mult_Alt", "")
        if "SnpEff_WARNING/ERROR" in df.columns.values:
            df.rename(columns={"SnpEff_WARNING/ERROR": "WARNING/ERROR"}, inplace=True)
        if "WARNING/ERROR" in df.columns.values:
            cols_to_group = self._rows + ["WARNING/ERROR"]
        else:
            cols_to_group = self._rows
        
        grouped = df.groupby(cols_to_group)
        
        df.set_index(cols_to_group, inplace=True)

        self.label_mult_alts(grouped, df)
        
        df.reset_index(inplace=True)
     
    def label_mult_alts(self, grouped, df):
        preliminary_dict = {}
        multalt_dict     = {}
        for k, gp in grouped:
            k_list = []
            for item in k:
                k_list.append(item)
            
            val = str(k_list[3]) ##CHROM,POS,REF,GENE_SYMBOL,[WARNING/ERROR]
            del k_list[3]
            key = "_".join(k_list) ##ANNOTATED_ALLELE
           
            if key in preliminary_dict:
                if key in multalt_dict:
                    multalt_dict[key].append(val)
                else:
                    multalt_dict[key] = [preliminary_dict[key], val]
            else:
                preliminary_dict[key] = val

        for key, vals in multalt_dict.iteritems():
            split_key = key.split("_")
            chrom = split_key[0]
            pos = split_key[1]
            ref = split_key[2]
            gene_symbol = split_key[3]
            warn_err = split_key[4]
            
            for val in vals:
                df.loc[(str(chrom), str(pos), str(ref), str(val), str(gene_symbol), str(warn_err)), "Mult_Alt"] = "True"

    def _is_compatible(self, dataframe):      
        unpivoted_df = self._transform(dataframe)

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
     
        
        
class EpeeVariantPivoter():
    COLUMNS = ["SAMPLE_NAME"]
 
    @staticmethod
    def _exclude_errors_and_warnings(df):
        try:
            filtered_df = df[df['WARNING/ERROR']=='.']
            return filtered_df
        except KeyError as e:
            raise PivotError("File is missing SnpEff_WARNING/ERROR column; review input files")
        else:
            raise
            
    @staticmethod
    def _build_transform_method(rows, columns, pivot_values):
        def transform(df):
            if "SnpEff_WARNING/ERROR" in df.columns.values:
                df.rename(columns={"SnpEff_WARNING/ERROR": "WARNING/ERROR"}, inplace=True)

            filtered_df = EpeeVariantPivoter._exclude_errors_and_warnings(df)

            expanded_df = expand_format(filtered_df, pivot_values, rows)

            return project_prepivot(expanded_df, pivot_values, rows, columns)
        return transform
            
    def __init__(self, input_keys, pivot_values):
        self._pivoter = VariantPivoter(
            input_keys, 
            EpeeVariantPivoter.COLUMNS, 
            pivot_values,
            EpeeVariantPivoter._build_transform_method(input_keys, EpeeVariantPivoter.COLUMNS, pivot_values))

    def pivot(self):
        return self._pivoter.pivot()
        
    def add_file(self, path, reader, header_index):
        return self._pivoter.add_file(path, reader, header_index)
        
    def is_compatible(self, initial_df):
        try:
            self._pivoter._is_compatible(initial_df)
            return True
        except PivotError:
            return False
       
class VcfVariantPivoter():
    COLUMNS = ["SAMPLE_NAME"]
   
    @staticmethod
    def _build_transform_method(rows, columns, pivot_values):
        def transform(df):
            expanded_df = expand_format(df, pivot_values, rows)
            
            return project_prepivot(expanded_df, pivot_values, rows, columns)
        return transform
        
    def __init__(self, input_keys, pivot_values):
        self._pivoter = VariantPivoter(
        input_keys, 
        VcfVariantPivoter.COLUMNS, 
        pivot_values,
        VcfVariantPivoter._build_transform_method(input_keys, VcfVariantPivoter.COLUMNS, pivot_values))
        
    def pivot(self):
        return self._pivoter.pivot()
    
    def add_file(self, path, reader, header_index):
        self._pivoter.add_file(path, reader, header_index)

    def is_compatible(self, initial_df):
        try:
            self._pivoter._is_compatible(initial_df)
            return True
        except PivotError:
            return False

def validate_parameters():
    invalid_fields = []
    
    fields = header_names.split("\t")
    for key in input_keys:
        if key not in fields:
            invalid_fields.append(key)
    
    invalid_tags = validate_format_tags(fields)
                
    message = "Invalid input parameter(s) "
    raise_err = 0
    if invalid_fields != []:
        message += str(invalid_fields)
        raise_err = 1
    if invalid_tags != []:
        message += str(invalid_tags)
        raise_err = 1
    
    if raise_err == 1:
        raise PivotError(message)

def validate_format_tags(fields):
    invalid_tags = []
    for line in first_line:
        my_line = line.split("\t")
        for index, field in enumerate(fields):
            if field == "FORMAT":
                format = my_line[index]
                format_tags = format.split(":")
                
                for val in pivot_values:
                    if val not in format_tags:
                        invalid_tags.append(val)
    
    return invalid_tags
    
def build_pivoter(path, reader, input_keys, pivot_values, header_index):
    pivoters = [EpeeVariantPivoter(input_keys, pivot_values), VcfVariantPivoter(input_keys, pivot_values)]
    initial_df  = create_initial_df(path, reader, header_index)

    for pivoter in pivoters:
        if pivoter.is_compatible(initial_df):
            return pivoter

    raise PivotError("Input file is not compatible with defined pivoters; review input file or contact system admin.")
    
def create_initial_df(path, reader, header_index):
    initial_df = pd.read_csv(reader, sep="\t", header=header_index, dtype='str')

    return initial_df
  
def append_to_annot_df(initial_df, annot_df):  
    format_index = initial_df.columns.get_loc("FORMAT")
    last_column = initial_df.shape[1]
    
    annotation_columns = []
    annotation_columns.extend(initial_df.columns)
    
    ##remove column if it is format column or any column to the right of format
    for i in range(format_index, last_column):
        annotation_columns.remove(initial_df.columns[i])
    if "INFO" in annotation_columns:
        annotation_columns.remove("INFO")
    
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
  
def melt_samples(df): 
    first_sample_column_index = df.columns.get_loc("FORMAT") + 1
    first_sample_column = df.ix[:, first_sample_column_index]
    last_column = df.shape[1]
    
    first_sample_name = df.columns[first_sample_column_index].split("_")

    ##find sample names -- a bit naive because it just looks to see if subsequent fields are structured the same as first sample field. this will likely break if sample name is not structured like "Sample_1"
    sample_names = []
    for i in range(first_sample_column_index, last_column):
        field_name = df.columns[i]
        split_field_name = field_name.split("_")
        
        if len(first_sample_name) == len(split_field_name):
            sample_names.append(field_name)

    column_list = []
    column_list.extend(df.columns)

    for col in sample_names:
        if col in column_list: 
            column_list.remove(col)
    
    melted_df = melt(df, id_vars=column_list, var_name='SAMPLE_NAME', value_name='SAMPLE_DATA')
    
    return melted_df
    
def expand_format(df, formats_to_expand, rows):     
    df = melt_samples(df)

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
    
    pd.pivot_table(joined_df, rows=rows+["SAMPLE_NAME"], cols="FORMAT2", values="VALUE2", aggfunc=lambda x: x)
    
    try:
        pivoted_df = pd.pivot_table(joined_df, rows=rows+["SAMPLE_NAME"], cols="FORMAT2", values="VALUE2", aggfunc=lambda x: x)
    except:
        # print "ugh"
        raise PivotError("Cannot pivot data.")

    for tag in formats_to_expand:
        if tag not in pivoted_df.columns:
            raise PivotError("Error: specified format tag {0} was not found in the file\n".format(tag))

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

def merge_samples(reduced_df, combined_df):
    if combined_df.empty:
        combined_df = reduced_df
    else:
        if len(combined_df.columns) == len(reduced_df.columns):
            combined_df = combined_df.append(reduced_df, ignore_index=True)
        else:
            raise PivotError("Columns do not match across files")

    return combined_df    
    
def rearrange_columns(output, pivot_values):
    ##change tuples to strings:
    lst = []
    for i in list(output.columns.values):
        if type(i) is tuple:
            i = "_".join(i)
        lst.append(i)
    output.columns = lst 
    
    ##change order of columns:
    lst = change_order(lst, pivot_values)    
    output = output.ix[:,lst]

    return output
    
def change_order(lst, pivot_values):
    all_pivot_values_lst = []
    for pivot_value in pivot_values:
        pivot_value_lst = []
        for i, header in enumerate(lst):    
            if re.search(pivot_value, header):
                pivot_value_lst.append(header)
        
        all_pivot_values_lst.extend(pivot_value_lst)
    
    meta_lst = []
    annot_lst = []
    for i, header in enumerate(lst):    
        if re.search("CHROM|POS|ID|REF|ALT|ANNOTATED_ALLELE|GENE_SYMBOL|Mult_Alt|IGV|UCSC", header):
            meta_lst.append(header)
        elif header not in all_pivot_values_lst:
            annot_lst.append(header) 
                     
    lst = meta_lst + all_pivot_values_lst + annot_lst

    return lst
    
def insert_links(joined_output):
    for row, col in joined_output.T.iteritems():
        joined_output.loc[row, "UCSC"] = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr" + joined_output.loc[row, "CHROM"] + ":" + str(int(joined_output.loc[row, "POS"]) - 50)+ "-" + str(int(joined_output.loc[row, "POS"]) + 50)
        
        joined_output.loc[row, "IGV"] = "localhost:60151/goto?locus=chr" + joined_output.loc[row, "CHROM"] + ":" + joined_output.loc[row, "POS"]
 
def style_workbook(output_path):
    wb = load_workbook(output_path)
    ws = wb.active
    
    for col in ws.columns:  
        if re.search("IGV", col[0].value):
            format_links(ws, col, "IGV")
        elif re.search("UCSC", col[0].value): 
            format_links(ws, col, "UCSC")
    
        count = 0
        colors = ["FFFF99", "FFCC00", "FF9900", "FF6600"]
        for tag in pivot_values:
            if re.search(tag, col[0].value):
                fill_cell(col, colors[count])
            count += 1
        
        if re.search("CHROM|POS|ID|REF|ALT|Mult_Alt|ANNOTATED_ALLELE|GENE_SYMBOL|IGV|UCSC", col[0].value):
            fill_cell(col, "C0C0C0")   

        if col[0].style.fill.start_color.index == "FFFFFFFF":
            fill_cell(col, "C9F2C9")
        # if re.search("Mult_Alt", col[0].value):
            # for row in col:
                # print row.value
    wb.save(output_path)
 
def format_links(ws, column, cell_value):
    for pos in column:
        if pos != column[0]:
            desired_cell = str(pos).strip("(<>)").split(".")[-1]
            ws[desired_cell].hyperlink = pos.value
            pos.value = cell_value
            pos.style.font.color.index = Color.BLUE
            pos.style.font.underline = Font. UNDERLINE_SINGLE
            
def fill_cell(column, color):
    column[0].style.fill.fill_type = Fill.FILL_SOLID
    column[0].style.fill.start_color.index = color
    
def process_files(sample_file_readers, pivot_builder=build_pivoter):
    first_file_reader = sample_file_readers[0]
    first_path        = str(first_file_reader)
    first_reader      = first_file_reader
    
    print "validating command line parameters"
    validate_parameters()
    
    print "determining file type"
    pivoter  = pivot_builder(first_path, first_reader, input_keys, pivot_values, header_index)
    # pivoter = EpeeVariantPivoter(pivot_values)
    
    for file_reader in sample_file_readers:
        path   = str(file_reader)
        reader = file_reader
        
        print "Reading: " + path
         
        pivoter.add_file(path, reader, header_index)
  
    print "validating annotations"
    pivoter._pivoter.validate_annotations()
    
    print "validating sample data"
    pivoter._pivoter.validate_sample_data()
  
    print "pivoting data"
    pivoted_df = pivoter.pivot()
    pivoted_df = pivoted_df.fillna("")
    
    print "joining sample data with annotations"
    joined_output = pivoter._pivoter.join_dataframes(pivoted_df)
    insert_links(joined_output)
  
    print "calculating mult alts"
    pivoter._pivoter.insert_mult_alt(joined_output)
    
    try:
        joined_output["CHROM"] = joined_output["CHROM"].apply(lambda x: int(x))
        joined_output["POS"] = joined_output["POS"].apply(lambda x: int(x))
    except:
        pass

    print "rearranging columns"
    complete_output = rearrange_columns(joined_output, pivot_values)
    
    print "sorting rows"
    sorted_df = pivoter._pivoter.sort_rows(complete_output)

    if "index" in sorted_df:
        del sorted_df["index"]

    print "writing to excel file: {0}".format(output_path)
    writer = ExcelWriter(output_path)
    sorted_df.to_excel(writer, "Variant_output", index=False, merge_cells = 0)  
    
    # print "writing to csv file: {0}".format(output_path)
    # sorted_df.to_csv(output_path, index=False, sep=",")  
    # print "almost done..."
   
    print "saving file"
    writer.save() 
    
    # print "styling workbook"
    # style_workbook(output_path)
   
    print "done"
    
if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pd.set_option('chained_assignment', None)
    
    if len(sys.argv) != 5:
        print "Invalid arguments."
        print "Usage: pivot.py [in_dir] [out_excel_file_path] [columns_for_pivot_keys] [format_tag(s)]"
        print "Example: pivot.py Rhim_input Rhim_output/output.xlsx CHROM,POS,REF,ALT GT,DP"
        exit(1)
    else:
        input_dir    = os.path.abspath(sys.argv[1])
        output_path  = os.path.abspath(sys.argv[2])
        input_keys   = sys.argv[3].split(",")
        pivot_values = sys.argv[4].split(",")

    output_dir, outfile_name = os.path.split(output_path)
    
    if not os.path.isdir(input_dir):
        print "Error. Specified input directory {0} does not exist".format(input_dir)
        exit(1)
        
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    sample_file_readers = []
    headers = []
    header_names = []
    first_line = []
    for item in listdir(input_dir):
        if isfile(join(input_dir, item)):
            f = open(input_dir + "/" + item, 'r')
            for num, line in enumerate(f):
                if line.startswith("#"):
                    headers.append(num)
                    header_names.append(line)
                else:
                    first_line.append(line)
                    break
            f.close()
            sample_file_readers.append(input_dir + "/" + item)

    header_index = headers[-1]
    header_names = header_names[-1]
    
    process_files(sample_file_readers)
    
    

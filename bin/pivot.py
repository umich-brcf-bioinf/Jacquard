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


class PivotError(Exception):
    """Base class for exceptions in this module."""
    pass 


class VariantPivoter():

    MISSING_REQUIRED_COLUMNS_ERROR="The columns of specified dataframe do not " +\
        "contain required columns {0}; review input data."
    NOOP_TRANSFORM = lambda x : x

    def __init__(self, rows, cols, pivot_values, transform=NOOP_TRANSFORM):
        self._rows = rows
        self._cols = cols
        self._pivot_values = pivot_values
        self._transform = transform
        self._combined_df = pd.DataFrame()

    def pivot(self):
        pivoted_df = pd.pivot_table(
            self._combined_df, values=self._pivot_values, rows=self._rows, cols=self._cols, 
            aggfunc=lambda x: x)

        try:
            pivoted_df = pivoted_df.applymap(lambda x: int(x))
        except:
            pass

        return pivoted_df
        
    def add_file(self, path, reader):
        initial_df  = create_initial_df(path, reader)
        unpivoted_df = self._is_compatible(initial_df)
        self._combined_df = merge_samples(unpivoted_df, self._combined_df)
        
    def join_dataframes(self, pivot_df, annot_df):
        output = merge(annot_df, pivot_df, left_on=self._rows, right_index=True, how='outer', sort=False)
        output = output.fillna('.')

        return output, annot_df
    
    def sort_rows(self, df):
        df.reset_index(inplace=True)
            
        try:
            df["CHROM"] = df["CHROM"].apply(lambda x: x.replace("chr", ""))
        except:
            pass
        
        df["CHROM"] = df["CHROM"].apply(lambda x: int(x))
        sorted_df = df.sort(self._rows)
        sorted_df = insert_mult_alt(sorted_df)
        
        grouped = sorted_df.groupby(self._rows)
        sorted_df.set_index(self._rows, inplace=True)
        sorted_df = self.label_mult_alts(grouped, sorted_df)

        return sorted_df
     
    def label_mult_alts(self, grouped, sorted_df):
        multalt_dict = {}
        for k, gp in grouped:
            key = str(k[0]) + "_" + str(k[1]) + "_" + str(k[2]) + "_" + str(k[4])
            val = str(k[3])
            if key in multalt_dict:
                multalt_dict[key].append(val)
            else:
                multalt_dict[key] = [val]
        
        for key, vals in multalt_dict.iteritems():
            if len(vals) != 1:
                split_key = key.split("_")
                chrom = int(split_key[0])
                pos = int(float(split_key[1]))
                ref = split_key[2]
                
                for val in vals:
                    for i, j in sorted_df.T.iteritems():
                        if chrom == i[0] and pos == i[1] and ref == i[2] and val == i[3]:
                            sorted_df.loc[i, "Mult_Alt"] = "True"
                            
        return sorted_df
    
    def _check_required_columns_present(self, dataframe):
        required_columns = set(self._rows + self._cols)

        if not required_columns.issubset(dataframe.columns.values):
            raise PivotError("Missing required columns; contact sysadmin.")
         
    def _check_pivot_is_unique(self, dataframe):
        group = self._rows + self._cols
        grouped_df = dataframe.groupby(group)
        
        if len(grouped_df.groups) != len(dataframe):
            raise PivotError("Duplicate keys would result in an invalid pivot; contact sysadmin.")
     
    def _is_compatible(self, dataframe):     
        unpivoted_df = self._transform(dataframe)
        
        self._check_required_columns_present(unpivoted_df)
        self._check_pivot_is_unique(unpivoted_df)
        
        return unpivoted_df
        
        
class EpeeVariantPivoter():
    ROWS = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
    COLUMNS = ["SAMPLE_NAME"]
 
    @staticmethod
    def _exclude_errors_and_warnings(df):
        try:
            # filtered_df = df.loc[df['WARNING/ERROR']=='.', ]
            # filtered_df = df[df['WARNING/ERROR']=='.']
            filtered_df = df[df["WARNING/ERROR"].isin(["."])]
            return filtered_df
        except KeyError as e:
            raise PivotError("File is missing WARNING/ERROR column; review input files")
        else:
            raise
            
    @staticmethod
    def _build_transform_method(rows, columns, pivot_values):
        def transform(df):
            filtered_df = EpeeVariantPivoter._exclude_errors_and_warnings(df)
            expanded_df = expand_format(filtered_df, "FORMAT", "Sample1", pivot_values)
            return project_prepivot(expanded_df, pivot_values, rows, columns)
        return transform
            
    def __init__(self, pivot_values):
        self._pivoter = VariantPivoter(
            EpeeVariantPivoter.ROWS, 
            EpeeVariantPivoter.COLUMNS, 
            pivot_values,
            EpeeVariantPivoter._build_transform_method(EpeeVariantPivoter.ROWS, EpeeVariantPivoter.COLUMNS, pivot_values))

    def pivot(self):
        return self._pivoter.pivot()
        
    def add_file(self, path, reader):
        return self._pivoter.add_file(path, reader)
        
    def is_compatible(self, initial_df):
        try:
            self._pivoter._is_compatible(initial_df)
            return True
        except PivotError:
            return False
       
class VcfVariantPivoter():
    ROWS = ["CHROM", "POS", "REF", "ALT"]
    COLUMNS = ["SAMPLE_NAME"]
    
    def __init__(self, pivot_values):
        self._pivoter = VariantPivoter(VcfVariantPivoter.ROWS, VcfVariantPivoter.COLUMNS, pivot_values)
        self._rows = VcfVariantPivoter.ROWS
        self._cols = VcfVariantPivoter.COLUMNS
        
    def pivot(self):
        return self._pivoter.pivot()

    def add_file(self, path, reader):
        self._pivoter.add_file(path, reader)

    def is_compatible(self, initial_df):
        try:
            self._pivoter._is_compatible(initial_df)
            return True
        except PivotError:
            return False


def build_pivoter(path, reader, pivot_values):
    pivoters = [EpeeVariantPivoter(pivot_values), VcfVariantPivoter(pivot_values)]
    initial_df  = create_initial_df(path, reader)
    
    for pivoter in pivoters:
        if pivoter.is_compatible(initial_df):
            return pivoter

    raise PivotError("Input file is not compatible with defined pivoters; review input file or contact system admin.")

        
def create_initial_df(path, reader):
    # initial_df = pd.read_csv(reader, sep="\t", header=False, dtype='str') ##for testing
    initial_df = pd.read_csv(reader, sep="\t", header=1, dtype='str') ##for actual epee data
    
    sample_name = os.path.basename(path)
    if "SAMPLE_NAME" not in initial_df:
        initial_df['SAMPLE_NAME'] = sample_name
    
    return initial_df

def expand_format(df, format_column_name, value_column_name, formats_to_expand):
    for row, col in df.T.iteritems():
        format_lst = col[format_column_name].split(":")
        sample_lst = col[value_column_name].split(":")
    
        format_dict = dict(zip(format_lst, sample_lst))
        for item in formats_to_expand:
            if item in format_dict:
                df.loc[row, item] = str(format_dict[item])
                
    return df
    
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
    
def rearrange_columns(output, sheet):
    ##change tuples to strings:
    lst = []
    for i in list(output.columns.values):
        if type(i) is tuple:
            i = "_".join(i)
        lst.append(i)
    output.columns = lst 
    
    ##change order of columns:
    lst = change_order(lst, sheet)    
    output = output.ix[:,lst]

    return output
    
def change_order(lst, sheet):
    if sheet == "variant":
        meta_lst = []
        freq_lst = []
        dp_lst = []
        annot_lst = []
        for i, header in enumerate(lst):    
            if re.search("CHROM|POS|REF|ANNOTATED_ALLELE|GENE_SYMBOL", header):
                meta_lst.append(header)

            elif re.search("FREQ", header):
                freq_lst.append(header)
            elif re.search("COVERAGE", header):
                dp_lst.append(header)
            else:
                annot_lst.append(header)    
        
        lst = meta_lst + freq_lst + dp_lst + annot_lst
        
    elif sheet == "gene":
        gene_lst = []
        level_lst = []
        impact_lst = []
        rank_lst = []
        for i, header in enumerate(lst):    
            if re.search("GENE_SYMBOL", header):
                gene_lst.append(header)
            elif re.search("IMPACT", header):
                level_lst.append(header)
            elif re.search("Impact_Damaging", header):
                impact_lst.append(header)
            elif re.search("Combined", header):
                rank_lst.append(header)    
        
        lst = gene_lst + level_lst + impact_lst + rank_lst
    
    return lst
    
def insert_mult_alt(df):
    df.insert(0, "Mult_Alt", "")

    return df
    
def process_files(sample_file_readers, pivot_builder=build_pivoter):
    pivot_values = ["GT"]
    
    first_file_reader = sample_file_readers[0]
    first_path        = str(first_file_reader)
    first_reader      = first_file_reader
    
    pivoter  = pivot_builder(first_path, first_reader, pivot_values)
    
    for file_reader in sample_file_readers:
        path   = str(file_reader)
        reader = file_reader
        
        pivoter.add_file(path, reader)
  
    pivoted_df = pivoter.pivot()
    pivoted_df = pivoted_df.fillna('.')
  
    annot_path = script_dir + '/../combined_annotation.txt'
    df_annot = pd.read_csv(annot_path, header=False, sep="\t", dtype='str')

    joined_output, df_annot = pivoter._pivoter.join_dataframes(pivoted_df, df_annot)
    complete_output = rearrange_columns(joined_output, "variant")
    
    sorted_df = pivoter._pivoter.sort_rows(complete_output)

    if "index" in sorted_df:
        del sorted_df["index"]
    
    writer = ExcelWriter(script_dir + '/../output/test_output.xlsx')
    sorted_df.to_excel(writer, "Variant_output", index=True, merge_cells = 0)  
    writer.save() 

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    sample_path = script_dir + '/../Rhim_input'
    sample_file_readers = []
    
    for item in listdir(sample_path):
        if isfile(join(sample_path, item)):
            sample_file_readers.append(sample_path + "/" + item)

    process_files(sample_file_readers)

    
    
##OLD CODE from here down:    
def process_file(df_annot):
    # sample_path = script_dir + '\input\combined_All_Genes.txt'
    
    # sample_path = script_dir + '\\Malek_input'
    sample_path = script_dir + '/../Rhim_input'
    epee = 1
    annotation = 1
    # format = "RNA_FREQ"
    format = ["RNA_FREQ", "RNA_COVERAGE"]
    fname = "Rhim_output"
    
    # sample_path = script_dir + '\\Malek_922_input'
    # epee = 0
    # annotation = 0
    # format = "FREQ"
    # fname = "Malek_922_all_vcfs_output"
    
    count = 1
    complete_df = pd.DataFrame()
    complete_df_gs = pd.DataFrame()

    if epee:
        header = 1
    else:
        header = False
    
    for item in listdir(sample_path):
        if isfile(join(sample_path, item)):
            print "Reading file " + str(count) + ": " + item
            
            df = pd.read_csv(sample_path + "/" + item, header=header, sep=r"\t")
            df2 = pd.read_csv(sample_path + "/" + item, header=header, sep=r"\t")
            
            print df
            exit(1)
            
            if "SAMPLE_NAME" not in df:
                df["SAMPLE_NAME"] = item
                df2["SAMPLE_NAME"] = item
            
            df_var = pivot(df, epee, format)
            if complete_df.empty:
                complete_df = df_var
            else:
                complete_df = complete_df.join(df_var, how="outer")
            
            if epee:
                df_gs = pivot2(df2)
                if complete_df_gs.empty:
                    complete_df_gs = df_gs
                else:
                    complete_df_gs = complete_df_gs.join(df_gs, how="outer")

        else:
            print "Skipping " + item
        
        count += 1
    
    complete_df_gs = calculate_impact_score(complete_df_gs)
    sorted_df_gs = rearrange_columns(complete_df_gs, "gene")

    if annotation:
        joined_output, df_annot = join_dataframes(complete_df, df_annot, epee)
        complete_output = rearrange_columns(joined_output, "variant")
        
        sorted_df = sort_rows(complete_output, epee)
        
        if "index" in sorted_df:
            del sorted_df["index"]
            
        write_files(fname, epee, sorted_df, sorted_df_gs, df_annot)
    
    else:
        sorted_df = sort_rows(complete_df, epee)
        write_files(fname, epee, sorted_df, sorted_df_gs, df_annot)
        
    print "done"        
        
def pivot(df, epee, format):
    # df = df[df['WARNING/ERROR']=='.']

    # df, values = expand_format(df, "FORMAT", "Sample1", format)
    values = format
    if epee:
        df = pd.pivot_table(df, values=values, rows=["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"], cols=["SAMPLE_NAME"], aggfunc=lambda x: x)
    else:
        df = pd.pivot_table(df, values=values, rows=["CHROM", "POS", "REF", "ALT"], cols=["SAMPLE_NAME"], aggfunc=lambda x: x)    
    
    try:
        df = df.applymap(lambda x: int(x))
    except:
        pass
    
    return df
    
def pivot2(df2):        
    # if "HIGHEST_IMPACT" in df2:
        # df2 = df2[df2['HIGHEST_IMPACT'] != 'MODIFIER']
    
    df_level = pd.pivot_table(df2, values=["ID"], rows=["GENE_SYMBOL", "SAMPLE_NAME"], cols=["HIGHEST_IMPACT"], aggfunc=lambda x: len(x))

    df_level = df_level["ID"]
    df_level = df_level.fillna(0)

    df_level["IMPACT"] = ""
    
    if "HIGH" not in df_level:
        df_level["HIGH"] = 0
    if "MODERATE" not in df_level:
        df_level["MODERATE"] = 0
    if "LOW" not in df_level:
        df_level["LOW"] = 0
    if "MODIFIER" not in df_level:
        df_level["MODIFIER"] = 0
    
    for col, data in df_level.T.iteritems():
        df_level["IMPACT"][col] = str(int(data["HIGH"])) + "|" + str(int(data["MODERATE"])) + "|" + str(int(data["LOW"])) + "|" + str(int(data["MODIFIER"]))
        
    df_level["IMPACT_SCORE"] = 100000.0 *df_level["HIGH"] + df_level["MODERATE"] + df_level["LOW"]/100000.0 + df_level["MODIFIER"]/10**12
    
    del df_level["HIGH"]
    del df_level["LOW"]
    del df_level["MODERATE"]
    del df_level["MODIFIER"]
    
    df_level = df_level.unstack("SAMPLE_NAME")
    
    df_impact = pd.pivot_table(df2, values=["Impact_Damaging"], rows=["GENE_SYMBOL"], cols=["SAMPLE_NAME"], aggfunc=lambda x: sum(x))
    
    df_combined = df_level.join(df_impact, how='outer')
    
    return df_combined    

def calculate_impact_score(df):
    df["IMPACT_Rank"]            = df["IMPACT_SCORE"].sum(axis=1)
    df["Impact_Damaging_Rank"] = df["Impact_Damaging"].sum(axis=1)
    # df["Combined_Rank"]        = df["Impact_Damaging_Rank"].add(df["IMPACT_Rank"])
    
    df["Impact_Damaging"] = df["Impact_Damaging"].applymap(lambda x: "." if x == 0.0 else x)
    
    df = df.sort("IMPACT_Rank", ascending=0)
    df = df.sort("Impact_Damaging_Rank", ascending=0)
    # df = df.sort("Combined_Rank", ascending=0)
    
    df["IMPACT_Rank"]            = df["IMPACT_Rank"].rank(ascending=0, method="min")
    df["Impact_Damaging_Rank"] = df["Impact_Damaging_Rank"].rank(ascending=0, method="min")
    # df["Combined_Rank"]        = df["Combined_Rank"].rank(ascending=0, method="min")

    del df["IMPACT_SCORE"]
    
    return df
    
def join_dataframes_old(df1, df2, epee):
    if epee:
        output = merge(df2, df1, left_on=["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"], right_index=True, how='right', sort=False)
    else:
        output = merge(df2, df1, left_on=["CHROM", "POS", "REF", "ALT"], right_index=True, how='right', sort=False)
    
    output = output.fillna('.')
    
    return output, df2
    
def rearrange_columns_old(output, sheet):
    ##change tuples to strings:
    lst = []
    for i in list(output.columns.values):
        if type(i) is tuple:
            i = "_".join(i)
        lst.append(i)
    output.columns = lst 
    
    ##change order of columns:
    lst = change_order(lst, sheet)    
    output = output.ix[:,lst]

    return output

def change_order_old(lst, sheet):
    if sheet == "variant":
        meta_lst = []
        freq_lst = []
        dp_lst = []
        annot_lst = []
        for i, header in enumerate(lst):    
            if re.search("CHROM|POS|REF|ANNOTATED_ALLELE|GENE_SYMBOL", header):
                meta_lst.append(header)

            elif re.search("FREQ", header):
                freq_lst.append(header)
            elif re.search("COVERAGE", header):
                dp_lst.append(header)
            else:
                annot_lst.append(header)    
        
        lst = meta_lst + freq_lst + dp_lst + annot_lst
        
    elif sheet == "gene":
        gene_lst = []
        level_lst = []
        impact_lst = []
        rank_lst = []
        for i, header in enumerate(lst):    
            if re.search("GENE_SYMBOL", header):
                gene_lst.append(header)
            elif re.search("IMPACT", header):
                level_lst.append(header)
            elif re.search("Impact_Damaging", header):
                impact_lst.append(header)
            elif re.search("Combined", header):
                rank_lst.append(header)    
        
        lst = gene_lst + level_lst + impact_lst + rank_lst
    
    return lst
    
def sort_rows_old(df, epee):
    df.reset_index(inplace=True)
        
    try:
        df["CHROM"] = df["CHROM"].apply(lambda x: x.replace("chr", ""))
    except:
        pass
    
    df["CHROM"] = df["CHROM"].apply(lambda x: int(x))

    if epee:
        sorted_df = df.sort(["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"])
    else:
        sorted_df = df.sort(["CHROM", "POS", "REF", "ALT"])
        
    # sorted_df["CHROM"] = sorted_df["CHROM"].apply(lambda x: "chr" + str(x))

    # pos_array = []
    # mult_alts = []
    # for pos in sorted_df["POS"]:
        # if pos in pos_array:
            # mult_alts.append(pos)
        # else:
            # pos_array.append(pos)
            
    sorted_df = insert_mult_alt(sorted_df)

    ### for position in mult_alts:
        ### sorted_df.loc[sorted_df["POS"] == position, "Mult_Alt"] = "True"

    grouped = sorted_df.groupby(["CHROM", "POS", "REF", "ANNOTATED_ALLELE"])
    
    if epee:
        sorted_df.set_index(["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"], inplace=True)
    else:
        sorted_df.set_index(["CHROM", "POS", "REF", "ALT"], inplace=True)
    
    #label mult-alts
    multalt_dict = {}
    for k, gp in grouped:
        key = str(k[0]) + "_" + str(k[1]) + "_" + str(k[2])
        val = str(k[3])
        if key in multalt_dict:
            multalt_dict[key].append(val)
        else:
            multalt_dict[key] = [val]
    
    for key, vals in multalt_dict.iteritems():
        if len(vals) != 1:
            split_key = key.split("_")
            chrom = int(split_key[0])
            pos = int(split_key[1])
            ref = split_key[2]
            
            for val in vals:
                for i, j in sorted_df.T.iteritems():
                    if chrom == i[0] and pos == i[1] and ref == i[2] and val == i[3]:
                        sorted_df.loc[i, "Mult_Alt"] = "True"
        
    return sorted_df    
    
def insert_mult_alt_old(df):
    df.insert(0, "Mult_Alt", "")

    return df
    
def write_files(fname, epee, output, df_gs, df_annot):
    df_gs = df_gs.fillna('.')

    # writer = ExcelWriter(script_dir + '/output/output.xlsx')
    writer = ExcelWriter(script_dir + '/../output/' + fname + '.xlsx')
    
    output.to_excel(writer, "Variant_output", merge_cells = 0)    
    
    if epee:
        df_gs.to_excel(writer, "Gene_output", merge_cells = 0)    
        # df_annot.to_excel(writer, "Original", merge_cells = 0)    
        
    writer.save()    
    
# if __name__ == "__main__":
    # script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # annot_path = script_dir + '/../combined_annotation.txt'
    # df_annot = pd.read_csv(annot_path, header=False, sep=r"\t")
    # process_file(df_annot)
    
    # sample_path = script_dir + '\\..\\Rhim_input'
    # sample_file_readers = {}
    # sample_file_readers = []
    
    # for item in listdir(sample_path):
        # if isfile(join(sample_path, item)):
            # sample_file_readers[str(sample_path + "/" + item)] = sample_path + "/" + item
            # sample_file_readers.append(sample_path + "\\" + item)

    # process_files(sample_file_readers)
    
        
def pivot_samples():
    sample_path = script_dir + '\..\input\combined_All_Genes.txt'
    df = pd.read_csv(sample_path, header=False, sep=r"\t")
    # df = df[df['WARNING/ERROR']=='.']
    
    df = pd.pivot_table(df, values=["FORMAT", "Sample1"], rows=["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL", "WARNING/ERROR"], cols=["SAMPLE_NAME"], aggfunc=lambda x: x.iloc[0])  ##aggfunc=lambda x: x.iloc[0] takes the first item aggfunc=lambda x: x.iloc[0].split(':')

    for col, data in df.iteritems():
        sample_name = col[1]
        
        format_col = df["FORMAT"][sample_name].dropna()
        format_col = format_col.apply(lambda x: x.split(':'))
        format_tags = format_col[-1]
        
        sample_col = df["Sample1"][sample_name].dropna()
        sample_col = sample_col.apply(lambda x: x.split(':'))
        sample_data = sample_col[-1]
        
        format_dict = dict(zip(format_tags, sample_data))
        
        # desired_tags = format_tags
        desired_tags = ["GT", "DP", "AD"]
        desired_format = []
        desired_data = []
        
        for desired_tag in desired_tags:
            if desired_tag in format_dict:
                desired_format.append(desired_tag)
                desired_data.append(format_dict[desired_tag])
        
        desired_format = ":".join(desired_format)
        desired_data = ":".join(desired_data)

        format_col = format_col.apply(lambda x: desired_format)
        df["FORMAT", sample_name] = format_col
        
        sample_col = sample_col.apply(lambda x: desired_data)
        df["Sample1", sample_name] = sample_col
    
    del df["FORMAT"]
    
    return df    
    

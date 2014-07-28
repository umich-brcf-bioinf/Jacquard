#!/usr/bin/python2.7
from collections import OrderedDict
import glob
import pandas as pd
from pandas import *
import unittest
from pandas.util.testing import assert_frame_equal
import pandas.util.testing as tm
from StringIO import StringIO
import pprint
import os
from os import listdir
from os.path import isfile, join
from bin.merge import PivotError, VariantPivoter, merge_samples, create_initial_df, build_pivoter, validate_parameters, rearrange_columns, determine_input_keys, get_headers_and_readers, create_dict, cleanup_df, combine_format_columns, remove_non_jq_tags, add_all_tags, sort_format_tags, determine_merge_execution_context, print_new_execution_context

def dataframe(input_data, sep="\t", index_col=None):
    def tupelizer(thing):
        if isinstance(thing, str) and thing.startswith("(") and thing.endswith(")"):
            return ast.literal_eval(thing)
        return thing

    df = pd.read_csv(StringIO(input_data), sep=sep, header=False, dtype='str', index_col=index_col)
    new_cols = [tupelizer(col) for col in list(df.columns.values)]
    df.columns = pd.core.index.Index(new_cols)
    
    return df

class MergeTestCase(unittest.TestCase):
    def test_addFiles(self):
        rows = ['COORDINATE']

        pivoter = VariantPivoter(rows)
        sample_A_file = \
'''COORDINATE\tFORMAT\tSamp1
1\tDP:ESAF\t1:0.2
2\tDP:ESAF\t12:0.2
3\tDP:ESAF\t31:0.2
4\tDP:ESAF\t6:0.2'''
        sample_B_file = \
'''COORDINATE\tFORMAT\tSamp2
1\tDP:ESAF\t5:0.2
2\tDP:ESAF\t2:0.2
3\tDP:ESAF\t74:0.2
4\tDP:ESAF\t25:0.2'''

        pivoter.add_file(StringIO(sample_A_file), 0, "file1")
        pivoter.add_file(StringIO(sample_B_file), 0, "file2")
        
        actual_df = pivoter._combined_df
        actual_df.columns.names = [""]
        expected_string = \
'''COORDINATE\tfile1|FORMAT\tfile1|Samp1\tfile2|FORMAT\tfile2|Samp2
1\tDP:ESAF\t1:0.2\tDP:ESAF\t5:0.2
2\tDP:ESAF\t12:0.2\tDP:ESAF\t2:0.2
3\tDP:ESAF\t31:0.2\tDP:ESAF\t74:0.2
4\tDP:ESAF\t6:0.2\tDP:ESAF\t25:0.2'''
        expected_df = dataframe(expected_string)
        expected_df.columns.names = [""]

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_isCompatible_raiseIfMissingRequiredColumns(self): 
        rows = ['COORDINATE', 'foo']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        pivoter = VariantPivoter(rows)
        
        input_string = \
'''COORDINATE\tFORMAT\tsample_A\tsample_B
1\tDP:ESAF\t10:0.2\t100:0.2
2\tDP:ESAF\t20:0.2\t200:0.2
3\tDP:ESAF\t30:0.2\t300:0.2
4\tDP:ESAF\t40:0.2\t400:0.2'''
        df = dataframe(input_string)
        
        self.assertRaises(PivotError, pivoter.is_compatible, df)
        
    def test_checkRequiredColumnsPresent(self):
        rows = ['Foo']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        pivoter = VariantPivoter(rows)
       
        expected_string = \
'''COORDINATE\tsample_A\tsample_B
1\t10\t100
2\t20\t200
3\t30\t300
4\t40\t400'''
        df = dataframe(expected_string)
        
        self.assertRaises(PivotError, pivoter._check_required_columns_present, df)
        
    def test_checkPivotIsUnique_DuplicateRaisesError(self):
        rows = ["Foo", "Bar", "Baz"]
        cols = ["Blah"]
        pivot_values = ['DP']
        pivoter = VariantPivoter(rows)
       
        expected_string = \
'''Foo\tBar\tBaz\tBlah
1\tA\t42\t2
1\tA\t42\t2'''

        df = dataframe(expected_string)

        self.assertRaises(PivotError, pivoter._check_pivot_is_unique, df)
        
    def test_validateSampleData_nonUniqueRows(self):
        rows = ["CHROM", "POS", "REF", "ALT"]

        input_string = \
'''CHROM\tPOS\tREF\tALT\tFORMAT\tSAMPLE_DATA
1\t2\tA\tCG\tGT\t10
1\t2\tA\tT\tGT\t13
1\t2\tA\tT\tGT\t12
2\t3\tC\tG\tGT\t11'''
        df = dataframe(input_string)
        combined_df = df
        
        pivoter = VariantPivoter(rows, combined_df)
        actual_df = pivoter.validate_sample_data()
        
        expected_string = \
'''CHROM\tPOS\tREF\tALT\tFORMAT\tSAMPLE_DATA
1\t2\tA\tCG\tGT\t10
1\t2\tA\tT\tGT\t^
2\t3\tC\tG\tGT\t11'''
        expected_df = dataframe(expected_string)

        tm.assert_frame_equal(expected_df, actual_df)
    
class PivotTestCase(unittest.TestCase):
    ##validate parameters
    def test_validateParameters_allValid(self):
        input_keys = ["CHROM", "POS", "REF"]
        first_line = ["1\t24\tA\tC\tGT:DP\tfoo;bar\t1/1:258"]
        header_names = "CHROM\tPOS\tREF\tALT\tFORMAT\tINFO\tsample2"

        output, message = validate_parameters(input_keys, first_line, header_names)

        self.assertEqual(0, output)
    
    def test_validateParameters_invalidKeys(self):
        input_keys = ["foo"]
        first_line = ["1\t24\tA\tC\tGT:DP\tfoo;bar\t1/1:258"]
        header_names = "CHROM\tPOS\tREF\tALT\tFORMAT\tINFO\tsample2"

        output, message = validate_parameters(input_keys, first_line, header_names)
        
        self.assertEqual(1, output)
        self.assertEqual("Invalid input parameter(s) ['foo']", message)
        
    def test_determineInputKeysVcf(self):
        input_dir = "test/test_input/test_input_keys_vcf"
        actual_lst = determine_input_keys(input_dir)
        
        expected_lst = ["CHROM", "POS", "REF", "ALT"]
        
        self.assertEquals(expected_lst, actual_lst)
        
    def test_determineInputKeysInvalid(self):
        input_dir = "test/test_input/test_input_keys_invalid"
        
        self.assertRaises(PivotError, determine_input_keys, input_dir)
    
    ##get headers, readers
    def test_getHeadersAndReaders(self):
        input_dir = "test/test_input/test_input_keys_txt"
        in_files = sorted(glob.glob(os.path.join(input_dir,"*")))
        sample_file_readers, headers, header_names, first_line, meta_headers = get_headers_and_readers(in_files)
        
        self.assertEquals([os.path.join(input_dir, "foo1.txt"), os.path.join(input_dir, "foo2.txt")], sample_file_readers)
        self.assertEquals([3,2], headers)
        self.assertEquals("CHROM\tPOS\tREF\tALT\tGENE_SYMBOL\tFORMAT\tSample_2384\tSample_2385\n", header_names)
        self.assertEquals(["1\t2342\tA\tT\tEGFR\tGT:DP\t1/1:241\t0/1:70\n", "1\t134\tG\tC\tEGFR\tGT:DP\t1/1:242\t0/1:546\n"], first_line)
        self.assertEquals(['##FORMAT=<ID=JQ_FOO'], meta_headers)
    def test_buildPivoter_invalidHeaderRaisesPivotError(self):
        input_string = \
'''COORDINATE\tFORMAT\tsample_A\tsample_B
1\tGT:ESAF\t10:0.2\t100:0.2
2\tGT:ESAF\t20:0.2\t200:0.2
3\tGT:ESAF\t30:0.2\t300:0.2
4\tGT:ESAF\t40:0.2\t400:0.2'''
        input_keys = ['CHROM', 'POS']
            
        self.assertRaises(PivotError, build_pivoter, StringIO(input_string), input_keys, 0)

    ##create_initial_df
    def test_createInitialDf(self):
        reader = StringIO( 
'''#CHROM\tPOS\tREF
1\t42\tA
2\t43\tC''');

        actual_df = create_initial_df(reader, 0)

        expected_data = StringIO(
'''#CHROM\tPOS\tREF
1\t42\tA
2\t43\tC''');

        expected_df = pd.read_csv(expected_data, sep="\t", header=False, dtype='str')

        tm.assert_frame_equal(expected_df, actual_df)

    def test_mergeSamples_emptyCombinedDf(self): 
        dataString = \
'''COORDINATE\tFORMAT\tsample_A\tsample_B
1\tGT:ESAF\t10:0.2\t100:0.2
2\tGT:ESAF\t20:0.2\t200:0.2
3\tGT:ESAF\t30:0.2\t300:0.2
4\tGT:ESAF\t40:0.2\t400:0.2'''
        df = pd.read_csv(StringIO(dataString), sep="\t", header=False, dtype='str')

        combined_df = pd.DataFrame()
        actual_df = merge_samples(df, combined_df, ["COORDINATE"])

        expected_data = StringIO(
'''COORDINATE\tFORMAT\tsample_A\tsample_B
1\tGT:ESAF\t10:0.2\t100:0.2
2\tGT:ESAF\t20:0.2\t200:0.2
3\tGT:ESAF\t30:0.2\t300:0.2
4\tGT:ESAF\t40:0.2\t400:0.2''')
        expected_df = pd.read_csv(expected_data, sep="\t", header=False, dtype='str')
        
        tm.assert_frame_equal(expected_df, actual_df)
    
    def test_mergeSamples_populatedCombinedDf(self): 
        dataString1 = \
'''COORDINATE\tfile1|FORMAT\tfile1|sample_A\tfile1|sample_B
1\tGT:ESAF\t10:0.2\t100:0.2
2\tGT:ESAF\t20:0.2\t200:0.2
3\tGT:ESAF\t30:0.2\t300:0.2
4\tGT:ESAF\t40:0.2\t400:0.2'''
        df1 = pd.read_csv(StringIO(dataString1), sep="\t", header=False, dtype='str')
        
        dataString2 = \
'''COORDINATE\tfile2|FORMAT\tfile2|sample_C\tfile2|sample_D
1\tGT:ESAF\t10:0.2\t100:0.2
2\tGT:ESAF\t20:0.2\t200:0.2
3\tGT:ESAF\t30:0.2\t300:0.2
4\tGT:ESAF\t40:0.2\t400:0.2'''
        df2 = pd.read_csv(StringIO(dataString2), sep="\t", header=False, dtype='str')

        combined_df = df2
        actual_df = merge_samples(combined_df, df1, ["COORDINATE"])

        expected_data = StringIO(
'''COORDINATE\tfile1|FORMAT\tfile1|sample_A\tfile1|sample_B\tfile2|FORMAT\tfile2|sample_C\tfile2|sample_D
1\tGT:ESAF\t10:0.2\t100:0.2\tGT:ESAF\t10:0.2\t100:0.2
2\tGT:ESAF\t20:0.2\t200:0.2\tGT:ESAF\t20:0.2\t200:0.2
3\tGT:ESAF\t30:0.2\t300:0.2\tGT:ESAF\t30:0.2\t300:0.2
4\tGT:ESAF\t40:0.2\t400:0.2\tGT:ESAF\t40:0.2\t400:0.2''')
        expected_df = pd.read_csv(expected_data, sep="\t", header=False, dtype='str')

        tm.assert_frame_equal(expected_df, actual_df)
    
    ##rearrange columns
    def test_rearrangeColumns(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tFORMAT\tsample_A\tsample_B\tINFO
1\t2\t.\tA\tG\tQUAL\tFILTER\tDP\t57\t57\tfoo
1\t3\t.\tC\tG\tQUAL\tFILTER\tDP\t58\t57\tfoo
8\t4\t.\tA\tT\tQUAL\tFILTER\tDP\t59\t57\tfoo
13\t5\t.\tT\tAA\tQUAL\tFILTER\tDP\t60\t57\tfoo
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tDP\t60\t57\tfoo
'''
        df = dataframe(input_string)
        actual_df = rearrange_columns(df)
        
        expected_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tDP\t57\t57
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tDP\t58\t57
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tDP\t59\t57
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tDP\t60\t57
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tDP\t60\t57
'''
        expected_df = dataframe(expected_string)

        tm.assert_frame_equal(expected_df, actual_df)
    
    
class CombineFormatTestCase(unittest.TestCase):
    def test_createDict(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tfile1|FORMAT\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tDP\t57\t57
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tDP\t58\t57
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tDP\t59\t57
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tDP\t60\t57
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tDP\t60\t57
'''
        df = dataframe(input_string)
        row = 0
        col = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "file1|FORMAT", "file1|sample_A", "file1|sample_B"]
        expected_file_dict = {'file1': [OrderedDict([('DP', '57'), ('sample_name', 'file1|sample_A')]), OrderedDict([('DP', '57'), ('sample_name', 'file1|sample_B')])]}
        file_dict, all_tags = create_dict(df, row, col)
        self.assertEqual(expected_file_dict, file_dict)
    
    def test_createDict_multFiles(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tfile1|FORMAT\tfile1|sample_A\tfile1|sample_B\tfile2|FORMAT\tfile2|sample_A\tfile2|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tDP\t57\t57\tDP\t57\t57
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tDP\t58\t57\tDP\t57\t57
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tDP\t59\t57\tDP\t57\t57
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tDP\t60\t57\tDP\t57\t57
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tDP\t60\t57\tDP\t57\t57
'''
        df = dataframe(input_string)
        row = 0
        col = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "file1|FORMAT", "file1|sample_A", "file1|sample_B", "file2|FORMAT", "file2|sample_A", "file2|sample_B"]
        expected_file_dict = {'file2': [OrderedDict([('DP', '57'), ('sample_name', 'file2|sample_A')]), OrderedDict([('DP', '57'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('DP', '57'), ('sample_name', 'file1|sample_A')]), OrderedDict([('DP', '57'), ('sample_name', 'file1|sample_B')])]}
        file_dict, all_tags = create_dict(df, row, col)
        self.assertEqual(expected_file_dict, file_dict)
        
#     def test_getConsensusFormatSample(self):
#         file_dict = {'file2': [OrderedDict([('DP', '57'), ('foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('DP', '57'), ('foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('DP', '57'), ('sample_name', 'file1|sample_A')]), OrderedDict([('DP', '57'), ('sample_name', 'file1|sample_B')])]}
#         all_tags = ["DP", "foo", "sample_name"]
#         
#         file_dict = get_consensus_format_sample(file_dict, all_tags)
#         expected_file_dict = {'file2': [OrderedDict([('DP', '57'), ('foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('DP', '57'), ('foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('DP', '57'), ('foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('DP', '57'), ('foo', '.'), ('sample_name', 'file1|sample_B')])]}
# 
#         self.assertEqual(expected_file_dict, file_dict)
        
    def test_cleanupDf(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile1|sample_A\tfile1|sample_B\tfile2|FORMAT\tfile2|sample_A\tfile2|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tDP:sample_name\t57:file1|sample_A\t57:file1|sample_B\tDP:sample_name\t57:file2|sample_A\t57:file2|sample_B
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tDP:sample_name\t58:file1|sample_A\t57:file1|sample_B\tDP:sample_name\t57:file2|sample_A\t57:file2|sample_B
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tDP:sample_name\t59:file1|sample_A\t57:file1|sample_B\tDP:sample_name\t57:file2|sample_A\t57:file2|sample_B
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tDP:sample_name\t60:file1|sample_A\t57:file1|sample_B\tDP:sample_name\t57:file2|sample_A\t57:file2|sample_B
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tnan:sample_name\tnan:file1|sample_A\tnan:file1|sample_B\tDP:sample_name\t57:file2|sample_A\t57:file2|sample_B
'''
        df = dataframe(input_string)
        file_dict = {"file1":"foo", "file2": "foo"}
        actual_df = cleanup_df(df, file_dict)
        
        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile1|sample_A\tfile1|sample_B\tfile2|sample_A\tfile2|sample_B
1\t2\t.\tA\tG\t.\t.\tfoo\tDP\t57\t57\t57\t57
1\t3\t.\tC\tG\t.\t.\tfoo\tDP\t58\t57\t57\t57
8\t4\t.\tA\tT\t.\t.\tfoo\tDP\t59\t57\t57\t57
13\t5\t.\tT\tAA\t.\t.\tfoo\tDP\t60\t57\t57\t57
13\t5\t.\tT\tAAAA\t.\t.\tfoo\t.\t.\t.\t57\t57
'''
        expected_df = dataframe(expected_string)
        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_cleanupDf_rearranged(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile1|sample_A\tfile1|sample_B\tfile2|FORMAT\tfile2|sample_A\tfile2|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tsample_name:DP\t57:file1|sample_A\t57:file1|sample_B\tDP:sample_name\t57:file2|sample_A\t57:file2|sample_B
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tsample_name:DP\t58:file1|sample_A\t57:file1|sample_B\tDP:sample_name\t57:file2|sample_A\t57:file2|sample_B
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tsample_name:DP\t59:file1|sample_A\t57:file1|sample_B\tDP:sample_name\t57:file2|sample_A\t57:file2|sample_B
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tsample_name:DP\t60:file1|sample_A\t57:file1|sample_B\tDP:sample_name\t57:file2|sample_A\t57:file2|sample_B
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tsample_name:nan\tnan:file1|sample_A\tnan:file1|sample_B\tDP:sample_name\t57:file2|sample_A\t57:file2|sample_B
'''
        df = dataframe(input_string)
        file_dict = {"file1":"foo", "file2": "foo"}
        actual_df = cleanup_df(df, file_dict)
        
        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile1|sample_A\tfile1|sample_B\tfile2|sample_A\tfile2|sample_B
1\t2\t.\tA\tG\t.\t.\tfoo\tDP\t57\t57\t57\t57
1\t3\t.\tC\tG\t.\t.\tfoo\tDP\t58\t57\t57\t57
8\t4\t.\tA\tT\t.\t.\tfoo\tDP\t59\t57\t57\t57
13\t5\t.\tT\tAA\t.\t.\tfoo\tDP\t60\t57\t57\t57
13\t5\t.\tT\tAAAA\t.\t.\tfoo\t.\t.\t.\t57\t57
'''
        expected_df = dataframe(expected_string)
        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_combineFormatColumns(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tfile1|FORMAT\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tJQ_DP\t57\t57
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tJQ_DP\t58\t57
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tJQ_DP\t59\t57
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tJQ_DP\t60\t57
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tnan\tnan\tnan
'''
        df = dataframe(input_string)
        actual_df = combine_format_columns(df)
        
        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tfile1|sample_A\tfile1|sample_B\tFORMAT
1\t2\t.\tA\tG\t.\t.\tfoo\t57\t57\tJQ_DP
1\t3\t.\tC\tG\t.\t.\tfoo\t58\t57\tJQ_DP
8\t4\t.\tA\tT\t.\t.\tfoo\t59\t57\tJQ_DP
13\t5\t.\tT\tAA\t.\t.\tfoo\t60\t57\tJQ_DP
13\t5\t.\tT\tAAAA\t.\t.\tfoo\t.\t.\t.
'''
        expected_df = dataframe(expected_string)
                
        print expected_df
        print actual_df
        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_combineFormatColumns_addTags(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tfile1|FORMAT\tfile1|sample_A\tfile1|sample_B\tfile2|FORMAT\tfile2|sample_A\tfile2|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t57:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t58:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t59:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t60:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tnan\tnan\tnan\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
'''
        df = dataframe(input_string)

        actual_df = combine_format_columns(df)
        
        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tfile1|sample_A\tfile1|sample_B\tfile2|sample_A\tfile2|sample_B\tFORMAT
1\t2\t.\tA\tG\t.\t.\tfoo\t0/1\t0/1\t0/1\t0/1\tJQ_GT
1\t3\t.\tC\tG\t.\t.\tfoo\t0/1\t0/1\t0/1\t0/1\tJQ_GT
8\t4\t.\tA\tT\t.\t.\tfoo\t0/1\t0/1\t0/1\t0/1\tJQ_GT
13\t5\t.\tT\tAA\t.\t.\tfoo\t0/1\t0/1\t0/1\t0/1\tJQ_GT
13\t5\t.\tT\tAAAA\t.\t.\tfoo\t.\t.\t0/1\t0/1\tJQ_GT
'''
        expected_df = dataframe(expected_string)
        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_removeNonJQTags(self):
        file_dict = {'file2': [OrderedDict([('DP', '57'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('DP', '57'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_B')])]}
        actual_file_dict = remove_non_jq_tags(file_dict)
        expected_file_dict = {'file2': [OrderedDict([('JQ_foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('JQ_foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('JQ_foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('JQ_foo', '.'), ('sample_name', 'file1|sample_B')])]}
       
        self.assertEquals(expected_file_dict, actual_file_dict)
    
    def test_removeNonJQTags_noTags(self):
        file_dict = {'file2': [OrderedDict([('DP', '57'), ('foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('DP', '57'), ('foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('DP', '57'), ('foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('DP', '57'), ('foo', '.'), ('sample_name', 'file1|sample_B')])]}
        actual_file_dict = remove_non_jq_tags(file_dict)
        expected_file_dict = {'file2': [OrderedDict([('sample_name', 'file2|sample_A')]), OrderedDict([('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('sample_name', 'file1|sample_A')]), OrderedDict([('sample_name', 'file1|sample_B')])]}
       
        self.assertEquals(expected_file_dict, actual_file_dict)
    
    def test_addAllTags(self):
        file_dict = {'file2': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('JQ_foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_B')])]}
        sample_keys = ["JQ_DP", "JQ_foo"]
        
        actual_file_dict = add_all_tags(file_dict, sample_keys)
        expected_file_dict = {'file2': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('JQ_foo', '1'), ('sample_name', 'file2|sample_B'), ('JQ_DP', '.')])], 'file1': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_B')])]}
        
        self.assertEquals(expected_file_dict, actual_file_dict)
        
    def test_sortFormatTags(self):
        file_dict = {'file2': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('JQ_foo', '1'), ('sample_name', 'file2|sample_B'), ('JQ_DP', '.')])], 'file1': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_B')])]}
        
        sorted_dict = sort_format_tags(file_dict)
        expected_sorted_dict = {'file2': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('JQ_DP', '.'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_B')])]}
        
        self.assertEquals(expected_sorted_dict, sorted_dict)
        
    def test_determineMergeContext(self):
        all_merge_context = []
        all_merge_column_context = []
        sample_columns = ["file1|samp1", "file1|samp2"]
        sample_file = "file1"
        count = 1
        all_merge_context, all_merge_column_context = determine_merge_execution_context(all_merge_context, all_merge_column_context, sample_columns, sample_file, count)
        
        self.assertEqual(["##jacquard.merge.file1=file1(['samp1', 'samp2'])"], all_merge_context)
        self.assertEqual(['##jacquard.merge.sample_column1=file1|samp1', '##jacquard.merge.sample_column2=file1|samp2'], all_merge_column_context)
        
    def test_printNewExecutionContext(self):
        out_file = MockWriter()
        execution_context = ["##jacquard", "##foo_bar", "##baz"]
        print_new_execution_context(execution_context, out_file)
        lines = out_file.lines()
        self.assertEqual(['##jacquard', '##foo_bar', '##baz'], lines)
        self.assertEqual(True, out_file.wasClosed)
        
class MockWriter():
    def __init__(self):
        self._content = []
        self.wasClosed = False

    def write(self, content):
        self._content.extend(content.splitlines())
        
    def lines(self):
        return self._content

    def close(self):
        self.wasClosed = True
        
#!/usr/bin/python2.7
import pandas as pd
from pandas import *
import unittest
from pandas.util.testing import assert_frame_equal
import pandas.util.testing as tm
from bin.pivot import PivotError, VariantPivoter, EpeeVariantPivoter, VcfVariantPivoter, pivot, expand_format, merge_samples, create_initial_df, project_prepivot, build_pivoter, append_to_annot_df, melt_samples, validate_parameters
from StringIO import StringIO
import pprint
import os
from os import listdir
from os.path import isfile, join

pd.set_option('chained_assignment', None)

def dataframe(input_data, sep="\t", index_col=None):
    return pd.read_csv(StringIO(input_data), sep=sep, header=False, dtype='str', index_col=index_col)
    
class VariantPivoterTestCase(unittest.TestCase):
    def test_validate_parameters_all_valid(self):
        input_keys = ["CHROM", "POS", "REF"]
        first_line = ["1\t24\tA\tC\tGT:DP\tfoo;bar\t1/1:258"]
        header_names = "CHROM\tPOS\tREF\tALT\tFORMAT\tINFO\tsample2"
        pivot_values = ["GT"]
        output, message = validate_parameters(input_keys, first_line, header_names, pivot_values)

        self.assertEqual(0, output)
    
    def test_validate_parameters_invalid_keys(self):
        input_keys = ["foo"]
        first_line = ["1\t24\tA\tC\tGT:DP\tfoo;bar\t1/1:258"]
        header_names = "CHROM\tPOS\tREF\tALT\tFORMAT\tINFO\tsample2"
        pivot_values = ["GT"]
        output, message = validate_parameters(input_keys, first_line, header_names, pivot_values)
        
        self.assertEqual(1, output)
        self.assertEqual("Invalid input parameter(s) ['foo']", message)
        
    def test_validate_parameters_invalid_pivot(self):
        input_keys = ["CHROM", "POS", "REF"]
        first_line = ["1\t24\tA\tC\tGT:DP\tfoo;bar\t1/1:258"]
        header_names = "CHROM\tPOS\tREF\tALT\tFORMAT\tINFO\tsample2"
        pivot_values = ["foo", "GT"]
        output, message = validate_parameters(input_keys, first_line, header_names, pivot_values)
        
        self.assertEqual(1, output)
        self.assertEqual("Invalid input parameter(s) ['foo']", message)

    def test_add_files(self):
        annot_df = pd.DataFrame()
        rows = ['COORDINATE']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        transform = lambda x : x
        pivoter = VariantPivoter(rows, cols, pivot_values, transform)
        sample_A_file = \
'''COORDINATE	FORMAT	Samp1	SAMPLE_NAME
1	DP	1	foo
2	DP	12	foo
3	DP	31	foo
4	DP	6	foo'''
        sample_B_file = \
'''COORDINATE	FORMAT	Samp1	SAMPLE_NAME
1	DP	5	foo
2	DP	2	foo
3	DP	74	foo
4	DP	25	foo'''

        pivoter.add_file("sample_A", StringIO(sample_A_file), 0)
        pivoter.add_file("sample_B", StringIO(sample_B_file), 0)
        actual_df = pivoter.pivot()

        actual_stringIO = StringIO()
        actual_df.to_csv(actual_stringIO, sep="\t")
        actual_flattened_df = pd.read_csv(StringIO(actual_stringIO.getvalue()), index_col=0, header=False, sep="\t", dtype='str')

        expected_string = \
'''COORDINATE	SAMPLE_NAME
1	foo
2	foo
3	foo
4	foo
1	foo
2	foo
3	foo
4	foo'''

        tm.assert_frame_equal(dataframe(expected_string), actual_flattened_df)
        
    def test_add_files_transformAsFilter(self):
        annot_df = pd.DataFrame()
        rows = ['COORDINATE']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        transform = lambda df : df[df["FILTER_COL"]=='.']
        pivoter = VariantPivoter(rows, cols, pivot_values, transform)
        sample_A_file = \
'''COORDINATE	FILTER_COL	SAMPLE_NAME	FORMAT
1	DP	.	foo
2	DP	.	foo
3	DP	Yep	foo
4	DP	.	foo
'''
        sample_B_file = \
'''COORDINATE	FILTER_COL	SAMPLE_NAME	FORMAT
1	.	foo	DP
2	Yep	foo	DP
3	.	foo	DP
4	.	foo	DP
'''
        pivoter.add_file("sample_A", StringIO(sample_A_file), 0)
        pivoter.add_file("sample_B", StringIO(sample_B_file), 0)
        actual_df = pivoter.pivot()

        actual_stringIO = StringIO()
        actual_df.to_csv(actual_stringIO, sep="\t")
        actual_flattened_df = pd.read_csv(StringIO(actual_stringIO.getvalue()), header=False, sep="\t", dtype='str')

        expected_string = \
'''	COORDINATE	SAMPLE_NAME
0	1	foo
2	3	foo
3	4	foo'''

        tm.assert_frame_equal(dataframe(expected_string), actual_flattened_df)

       
    def test_is_compatible_raiseIfMissingRequiredColumns(self): 
        rows = ['COORDINATE', 'foo']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        transform = lambda x : x
        pivoter = VariantPivoter(rows, cols, pivot_values, transform)
        
        input_string = \
'''COORDINATE	SAMPLE_NAME	sample_A	sample_B
1	foo	10	100
2	bar	20	200
3	foo	30	300
4	bar	40	400'''
        df = dataframe(input_string)
        
        self.assertRaises(PivotError, pivoter._is_compatible, df)
        
    def test_check_required_columns_present(self):
        rows = ['Foo']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        transform = lambda x : x
        pivoter = VariantPivoter(rows, cols, pivot_values, transform)
       
        expected_string = \
'''COORDINATE	sample_A	sample_B
1	10	100
2	20	200
3	30	300
4	40	400'''
        df = dataframe(expected_string)
        
        self.assertRaises(PivotError, pivoter._check_required_columns_present, df)
        
    def test_check_pivot_is_unique_DuplicateRaisesError(self):
        rows = ["Foo", "Bar", "Baz"]
        cols = ["Blah"]
        pivot_values = ['DP']
        transform = lambda x : x
        pivoter = VariantPivoter(rows, cols, pivot_values, transform)
       
        expected_string = \
'''Foo	Bar	Baz	Blah
1	A	42	2
1	A	42	2'''

        df = dataframe(expected_string)

        self.assertRaises(PivotError, pivoter._check_pivot_is_unique, df)
        
    def test_label_mult_alts(self):
        rows = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
        cols = ["SAMPLE_NAME"]
        pivot_values = ['GT']
        transform = lambda x : x
        pivoter = VariantPivoter(rows, cols, pivot_values, transform)
    
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
sample1	1	2	A	T	foo	.	GT	0/1	10	100
sample1	1	2	A	G	foo	.	GT	0/1	10	100
sample2	2	3	A	T	foo	.	GT	0/1	20	200'''
        df = dataframe(input_string)
        cols_to_group = rows + ["WARNING/ERROR"]
        grouped = df.groupby(cols_to_group)
        
        df.set_index(cols_to_group, inplace=True)
        pivoter.label_mult_alts(grouped, df)
        
        expected_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B	Mult_Alt
sample1	1	2	A	T	foo	.	GT	0/1	10	100	True
sample1	1	2	A	G	foo	.	GT	0/1	10	100	True
sample2	2	3	A	T	foo	.	GT	0/1	20	200	NaN'''
        expected_df = dataframe(expected_string)
        expected_df.set_index(cols_to_group, inplace=True)
        
        tm.assert_frame_equal(expected_df, df)
            
    def test_validate_sample_data_non_unique_cols(self):
        rows = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
        cols = ["SAMPLE_NAME"]
        pivot_values = ['GT']
        transform = lambda x : x
        
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	SAMPLE_NAME	SAMPLE_DATA
1	2	A	T	foo	.	GT	sampleA	13
1	2	A	T	foo	.	GT	sampleA	12
2	3	C	G	foo	.	GT	sampleB	11'''
        df = dataframe(input_string)
        combined_df = df
        
        pivoter = VariantPivoter(rows, cols, pivot_values, transform, combined_df)
        actual_df = pivoter.validate_sample_data()
        
        expected_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	SAMPLE_NAME	SAMPLE_DATA
1	2	A	T	foo	.	GT	sampleA	^
1	2	A	T	foo	.	GT	sampleA	^
2	3	C	G	foo	.	GT	sampleB	11'''
        expected_df = dataframe(expected_string)

        tm.assert_frame_equal(expected_df, actual_df)
        
    # not sure how to test if numpy.ndarray
    # def test_validate_sample_data_non_unique_rows(self):
        # rows = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
        # cols = ["SAMPLE_NAME"]
        # pivot_values = ['GT']
        # transform = lambda x : x
        
        # input_string = \
# '''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	SAMPLE_NAME	SAMPLE_DATA
# 1	2	A	T	foo	.	GT	sampleA	{0}
# 1	2	A	C	foo	.	GT	sampleA	11
# 2	3	C	G	foo	.	GT	sampleB	10'''.format(np.ndarray(shape=(1,2)))
        # df = dataframe(input_string)
        
        # pivoter = VariantPivoter(rows, cols, pivot_values, transform, df)
        # actual_df = pivoter.validate_sample_data()
        
        # expected_string = \
# '''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	SAMPLE_NAME	SAMPLE_DATA
# 1	2	A	T	foo	.	GT	sampleA	^
# 1	2	A	C	foo	.	GT	sampleA	11
# 2	3	C	G	foo	.	GT	sampleB	11'''
        # expected_df = dataframe(expected_string)
        # print actual_df.applymap(lambda x: type(x))
        
        # print actual_df
        # tm.assert_frame_equal(expected_df, actual_df)
    
    
    
class EpeeVariantPivoterTestCase(unittest.TestCase):
    def test_build_transform_method(self):
        input_df = dataframe(
            "#CHROM|POS|WARNING/ERROR|FORMAT|SAMPLE_1\n" + \
            "1|41|warn|GT:DP|0/41:410\n" + \
            "1|42|.|GT:DP|0/42:420\n",'|')
        transform = EpeeVariantPivoter._build_transform_method(['CHROM','POS'], ['SAMPLE_NAME'],['DP'])

        actual_df = transform(input_df)
        
        expected_df = dataframe(
            "CHROM|POS|SAMPLE_NAME|DP\n" + \
            "1|42|SAMPLE_1|420",sep='|')
            
        actual_df = actual_df.applymap(lambda x: str(x))
        
        tm.assert_frame_equal(expected_df, actual_df)


    def test_exclude_errors_warnings_throwsIfColumnMissing(self):
        input_df = dataframe(
            "#CHROM|POS\n" + \
            "1|A\n" + \
            "2|B",'|')
        input_keys = ['CHROM', 'POS', 'REF', 'ANNOTATED_ALLELE', 'GENE_SYMBOL']
        filter = EpeeVariantPivoter(input_keys, ["GT"])._pivoter._transform

        self.assertRaises(PivotError, filter, input_df)

    def test_is_compatible_trueIfCompatible(self):
        input_keys = ['CHROM', 'POS', 'REF', 'ANNOTATED_ALLELE', 'GENE_SYMBOL']
        pivot_values = ['GT']
        pivoter = EpeeVariantPivoter(input_keys, pivot_values)        
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
1	2	A	T	foo	.	GT	0/1	10	100
1	2	A	T	foo	WARN	GT	0/1	10	100
2	3	A	T	foo	.	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(True, pivoter.is_compatible(df))
    
    def test_is_compatible_falseIfMissingColumns(self):
        input_keys = ['CHROM', 'POS', 'REF', 'ANNOTATED_ALLELE', 'GENE_SYMBOL']
        pivot_values = ['DP']
        pivoter = EpeeVariantPivoter(input_keys, pivot_values)        
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
1	2	A	T	foo	.	GT	0/1	10	100
1	2	A	T	foo	WARN	GT	0/1	10	100
2	3	A	T	foo	.	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(False, pivoter.is_compatible(df))
        
    def test_is_compatible_trueIfWarningDuplicates(self):
        input_keys = ['CHROM', 'POS', 'REF', 'ANNOTATED_ALLELE', 'GENE_SYMBOL']
        pivot_values = ['GT']
        pivoter = EpeeVariantPivoter(input_keys, pivot_values)        
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
1	2	A	T	foo	.	GT	0/1	10	100
1	2	A	T	foo	WARN	GT	0/1	10	100
2	3	A	T	foo	.	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(True, pivoter.is_compatible(df))
        
    def test_is_compatible_falseIfDuplicates(self):
        input_keys = ['CHROM', 'POS', 'REF', 'ANNOTATED_ALLELE', 'GENE_SYMBOL']
        pivot_values = ['GT']
        pivoter = EpeeVariantPivoter(input_keys, pivot_values)        
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
1	2	A	T	foo	.	GT	0/1	10	100
1	2	A	T	foo	.	GT	0/1	10	100
2	3	A	T	foo	.	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(False, pivoter.is_compatible(df)) ##ERROR: function does not reduce

        
class VcfVariantPivoterTestCase(unittest.TestCase):
    def test_build_transform_method(self):
        input_df = dataframe(
            "#CHROM|POS|FORMAT|SAMPLE_1\n" + \
            "1|41|GT:DP|0/41:410\n" + \
            "1|42|GT:DP|0/42:420\n",'|')
        transform = VcfVariantPivoter._build_transform_method(['CHROM','POS'], ['SAMPLE_NAME'],['DP'])

        actual_df = transform(input_df)
        
        expected_df = dataframe(
            "CHROM|POS|SAMPLE_NAME|DP\n" + \
            "1|41|SAMPLE_1|410\n" + \
            "1|42|SAMPLE_1|420",sep='|')
            
        actual_df = actual_df.applymap(lambda x: str(x))

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_is_compatible_trueIfCompatible(self):
        input_keys = ['CHROM','POS','REF','ALT']
        pivot_values = ['GT']
        pivoter = VcfVariantPivoter(input_keys, pivot_values)        
        input_string = \
'''CHROM	POS	REF	ALT	GENE_SYMBOL	FORMAT	Sample1	sample_A	sample_B
1	2	A	T	foo	GT	0/1	10	100
2	3	A	T	foo	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(True, pivoter.is_compatible(df))  
        
    def test_is_compatible_falseIfDuplicates(self):
        input_keys = ['CHROM','POS','REF','ALT']
        pivot_values = ['DP']
        pivoter = VcfVariantPivoter(input_keys, pivot_values)        
        input_string = \
'''CHROM	POS	REF	ALT	GENE_SYMBOL	FORMAT	Sample1	sample_A	sample_B
1	2	A	T	foo	GT	0/1	10	100
1	2	A	T	foo	GT	0/1	10	100
2	3	A	T	foo	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(False, pivoter.is_compatible(df))
        
    def test_is_compatible_falseIfMissingColumns(self):
        input_keys = ['CHROM','POS','REF','ALT']
        pivot_values = ['DP']
        pivoter = VcfVariantPivoter(input_keys, pivot_values)        
        input_string = \
'''CHROM	POS	REF	ALT	GENE_SYMBOL	FORMAT	Sample1	sample_A	sample_B
1	2	A	T	foo	GT	0/1	10	100
2	3	A	T	foo	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(False, pivoter.is_compatible(df))
        
class PivotTestCase(unittest.TestCase):
    def test_build_pivoter_invalidHeaderRaisesPivotError(self):
        input_string = \
'''COORDINATE	FORMAT	sample_A	sample_B
1	GT	10	100
2	GT	20	200
3	GT	30	300
4	GT	40	400'''
        input_keys = ['CHROM', 'POS']
        pivot_values = ["GT"]
            
        self.assertRaises(PivotError, build_pivoter, "sampleA", StringIO(input_string), input_keys, pivot_values, 0)

    def test_build_pivoter_epeeHeader(self):
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
1	2	A	T	foo	WARN	GT	0/1	10	100
2	3	A	T	foo	.	GT	0/1	20	200
3	4	A	T	foo	.	GT	0/1	30	300
4	5	A	T	foo	.	GT	0/1	40	400'''
        input_keys = ['CHROM','POS','REF','ANNOTATED_ALLELE', 'GENE_SYMBOL']
        pivot_values = ["GT"]
        
        pivoter = build_pivoter("sampleA", StringIO(input_string), input_keys, pivot_values, 0)

        self.assertIsInstance(pivoter, EpeeVariantPivoter)
        
    def test_build_pivoter_vcfHeader(self):
        input_string = \
'''CHROM	POS	REF	ALT	FORMAT	SAMPLE_NAME	sample_A	sample_B
1	2	A	T	GT	foo	10	100
2	3	A	T	GT	foo	20	200
3	4	A	T	GT	bar	30	300
4	5	A	T	GT	bar	40	400'''
        input_keys = ['CHROM','POS','REF','ALT']
        pivot_values = ["GT"]
        
        pivoter = build_pivoter("sampleA", StringIO(input_string), input_keys, pivot_values, 0)
        
        self.assertIsInstance(pivoter, VcfVariantPivoter)
            
        
    ##create_initial_df
    def test_create_initial_df(self):
        reader = StringIO( 
'''#CHROM	POS	REF
1	42	A
2	43	C''');
        
        path = 'foo/bar/sampleA_AllGenes.txt'
        actual_df = create_initial_df(path, reader, 0)
        
        expected_data = StringIO(
'''#CHROM	POS	REF
1	42	A
2	43	C''');

        expected_df = pd.read_csv(expected_data, sep="\t", header=False, dtype='str')

        tm.assert_frame_equal(expected_df, actual_df)

        
    ##expand_format   
    def test_expand_format(self):
        dataString = \
'''CHROM	POS	FORMAT	Sample1
1	2	GT:DP	A:1
1	3	GT:DP	B:2
1	5	GT:DP	:3
13	3	GT:DP	D:'''
        input_df = pd.read_csv(StringIO(dataString), sep="\t", header=False)
        rows=["CHROM", "POS"]
        actual_df = expand_format(input_df, ["GT"], rows)
 
        expected_string = \
'''CHROM	POS	SAMPLE_NAME	DP	GT
1	2	Sample1	1	A
1	3	Sample1	2	B
1	5	Sample1	3	
13	3	Sample1		D'''

        expected_df = pd.read_csv(StringIO(expected_string), sep="\t", header=False, dtype={'DP':str})
        expected_df.fillna(value="", inplace=True)

        tm.assert_frame_equal(expected_df, actual_df, check_names=False)
        
    def test_expand_format_multipleFormats(self):
        dataString = \
'''CHROM	POS	FORMAT	Sample1
1	1	GT:DP	A:1
1	2	GT:DP	B:2
1	3	GT:DP	:3
1	4	GT:DP	D:'''
        input_df = pd.read_csv(StringIO(dataString), sep="\t", header=False)
        rows = ["CHROM", "POS"]
        actual_df = expand_format(input_df, ["GT", "DP"], rows)

        expectedString = \
'''CHROM\tPOS\tSAMPLE_NAME\tDP\tGT
1	1	Sample1	1	A
1	2	Sample1	2	B
1	3	Sample1	3	
1	4	Sample1		D'''
        expected_df = pd.read_csv(StringIO(expectedString), sep="\t", header=False, dtype={'DP':str})
        expected_df.fillna(value="", inplace=True)
     
        tm.assert_frame_equal(expected_df, actual_df, check_names=False)

    def test_expand_format_invalidFormatsIgnored(self):
        dataString = \
'''CHROM	POS	FORMAT	Sample1
1	1	GT:DP	A:1
1	2	GT:DP	B:2
1	3	GT:DP	:3
1	4	GT:DP	D:'''
        input_df = pd.read_csv(StringIO(dataString), sep="\t", header=False)
        rows = ["CHROM", "POS"]
      
        self.assertRaises(PivotError, expand_format, input_df, ["GT", "DP", "FOO"], rows)

        
    ##select_prepivot
    def test_project_prepivot(self):
        dataString = \
'''SAMPLE_NAME	CHROM	POS	REF	ALT	INFO	FORMAT	GT
sample1	chr1	1	A	T	blah	GT:DP	0/1
sample2	chr1	2	A	T	blah	GT:DP	0/1
sample3	chr1	3	A	T	blah	DP	1/1
sample6	chr1	4	A	T	blah	GT	1/1'''
        df = pd.read_csv(StringIO(dataString), sep="\t", header=False)
        pivot_values = ["GT"]
        rows = ["CHROM", "POS", "REF", "ALT"]
        columns = ["SAMPLE_NAME"]
        actual_df = project_prepivot(df, pivot_values, rows, columns)
        
        expected_dataString = \
'''SAMPLE_NAME	CHROM	POS	REF	ALT	GT
sample1	chr1	1	A	T	0/1
sample2	chr1	2	A	T	0/1
sample3	chr1	3	A	T	1/1
sample6	chr1	4	A	T	1/1'''
        expected_df = pd.read_csv(StringIO(expected_dataString), sep="\t", header=False)

        tm.assert_frame_equal(expected_df, actual_df)
        
    ##melt_samples
    def test_melt_samples(self):
        dataString = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	FORMAT	Sample_1	Sample_2
chr1	1	A	T	GT	0/1	1/1
chr1	2	A	T	GT	0/1	0/1
chr1	3	A	T	GT	1/1	1/1
chr1	4	A	T	GT	1/1	0/1'''
        df = pd.read_csv(StringIO(dataString), sep="\t", header=False, dtype='str')
        
        actual_df = melt_samples(df)
        
        expected_dataString = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	FORMAT	SAMPLE_NAME	SAMPLE_DATA
chr1	1	A	T	GT	Sample_1	0/1
chr1	2	A	T	GT	Sample_1	0/1
chr1	3	A	T	GT	Sample_1	1/1
chr1	4	A	T	GT	Sample_1	1/1
chr1	1	A	T	GT	Sample_2	1/1
chr1	2	A	T	GT	Sample_2	0/1
chr1	3	A	T	GT	Sample_2	1/1
chr1	4	A	T	GT	Sample_2	0/1'''
        expected_df = pd.read_csv(StringIO(expected_dataString), sep="\t", header=False, dtype='str')

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_melt_samples_trailing_field(self):
        dataString = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	FORMAT	Sample_1	Sample_2	Sample
chr1	1	A	T	GT	0/1	1/1	foo
chr1	2	A	T	GT	0/1	0/1	foo
chr1	3	A	T	GT	1/1	1/1	foo
chr1	4	A	T	GT	1/1	0/1	foo'''
        df = pd.read_csv(StringIO(dataString), sep="\t", header=False, dtype='str')
        
        actual_df = melt_samples(df)
        
        expected_dataString = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	FORMAT	Sample	SAMPLE_NAME	SAMPLE_DATA
chr1	1	A	T	GT	foo	Sample_1	0/1
chr1	2	A	T	GT	foo	Sample_1	0/1
chr1	3	A	T	GT	foo	Sample_1	1/1
chr1	4	A	T	GT	foo	Sample_1	1/1
chr1	1	A	T	GT	foo	Sample_2	1/1
chr1	2	A	T	GT	foo	Sample_2	0/1
chr1	3	A	T	GT	foo	Sample_2	1/1
chr1	4	A	T	GT	foo	Sample_2	0/1'''
        expected_df = pd.read_csv(StringIO(expected_dataString), sep="\t", header=False, dtype='str')

        tm.assert_frame_equal(expected_df, actual_df)
        
    ##merge_samples
    def test_merge_samples_emptyCombinedDf(self): 
        dataString = \
'''CHROM	POS	REF	ALT	SAMPLE_NAME	GT	DP
chr1	1	A	T	sample1	1/1	1
chr1	2	A	T	sample1	0/1	2
chr1	3	A	T	sample2	1/1	3
chr1	4	A	T	sample3	0/1	4'''
        df = pd.read_csv(StringIO(dataString), sep="\t", header=False, dtype='str')

        combined_df = pd.DataFrame()
        actual_df = merge_samples(df, combined_df)

        expected_data = StringIO(
'''CHROM	POS	REF	ALT	SAMPLE_NAME	GT	DP
chr1	1	A	T	sample1	1/1	1
chr1	2	A	T	sample1	0/1	2
chr1	3	A	T	sample2	1/1	3
chr1	4	A	T	sample3	0/1	4''')
        expected_df = pd.read_csv(expected_data, sep="\t", header=False, dtype='str')
        
        tm.assert_frame_equal(expected_df, actual_df)
    
    def test_merge_samples_populatedCombinedDf(self): 
        dataString1 = \
'''CHROM	POS	REF	ALT	SAMPLE_NAME	GT	DP
chr1	1	A	T	sample1	1/1	1
chr1	2	A	T	sample1	0/1	2
chr1	3	A	T	sample2	1/1	3
chr1	4	A	T	sample3	0/1	4'''
        df1 = pd.read_csv(StringIO(dataString1), sep="\t", header=False, dtype='str')
        
        dataString2 = \
'''CHROM	POS	REF	ALT	SAMPLE_NAME	GT	DP
chr1	1	A	T	sample4	0/1	5
chr1	2	A	T	sample5	0/1	2
chr2	3	A	T	sample6	1/1	3
chr3	4	A	T	sample7	0/1	4'''
        df2 = pd.read_csv(StringIO(dataString2), sep="\t", header=False, dtype='str')

        combined_df = df2
        
        actual_df = merge_samples(df1, combined_df)

        expected_data = StringIO(
'''CHROM	POS	REF	ALT	SAMPLE_NAME	GT	DP
chr1	1	A	T	sample4	0/1	5
chr1	2	A	T	sample5	0/1	2
chr2	3	A	T	sample6	1/1	3
chr3	4	A	T	sample7	0/1	4
chr1	1	A	T	sample1	1/1	1
chr1	2	A	T	sample1	0/1	2
chr1	3	A	T	sample2	1/1	3
chr1	4	A	T	sample3	0/1	4''')
        expected_df = pd.read_csv(expected_data, sep="\t", header=False, dtype='str')

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_append_to_annot_df(self):
        annot_df = pd.DataFrame()

        dataString = \
'''CHROM	POS	REF	ALT	FORMAT
chr1	1	A	T	DP:GT
chr1	2	A	T	DP:GT
chr1	3	A	T	DP:GT
chr1	4	A	T	DP:GT'''
        df = pd.read_csv(StringIO(dataString), sep="\t", header=False, dtype='str')
        
        annot_df = append_to_annot_df(df, annot_df)
        
        expected_dataString = \
'''CHROM	POS	REF	ALT
chr1	1	A	T
chr1	2	A	T
chr1	3	A	T
chr1	4	A	T'''
        expected_df = pd.read_csv(StringIO(expected_dataString), sep="\t", header=False, dtype='str')
        
        tm.assert_frame_equal(expected_df, annot_df)
        
    def test_annot_files_equal(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        annot_df = pd.DataFrame()
        
        file_list = [script_dir + "/test_input/P2_test_input.txt", script_dir + "/test_input/P5_test_input.txt"]
        for file in file_list:
            df = pd.read_csv(file, sep="\t", header=1, dtype='str', index_col=False)
            annot_df = append_to_annot_df(df, annot_df)

        expected_df = pd.read_csv(script_dir + "/test_annotation/P2P5_combined_annotation.txt", sep="\t", header=False, dtype='str', index_col=False)

        try:
            expected_df["CHROM"] = expected_df["CHROM"].apply(lambda x: x.replace("chr", ""))
            annot_df_df["CHROM"] = annot_df_df["CHROM"].apply(lambda x: x.replace("chr", ""))
           
            expected_df["CHROM"] = expected_df["CHROM"].apply(lambda x: int(x))
            annot_df_df["CHROM"] = annot_df_df["CHROM"].apply(lambda x: int(x))
        except:
            pass

        sort_by = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL", "WARNING/ERROR"]
        sorted_expected_df = expected_df.sort(sort_by)
        sorted_annot_df = annot_df.sort(sort_by)
        
        del sorted_expected_df['FORMAT']
        del sorted_expected_df['Sample1']
        del sorted_expected_df['P2_N_GT']
        del sorted_expected_df['P2_N_ESAF']
        del sorted_expected_df['INFO']

        sorted_expected_df.reset_index(inplace=True)
        sorted_annot_df.reset_index(inplace=True)
        
        if "index" in sorted_expected_df:
            del sorted_expected_df["index"]
       
        if "index" in sorted_annot_df:
            del sorted_annot_df["index"]

        sorted_annot_df.to_csv(script_dir + "/test_output/sorted_annot.txt", sep="\t")
        sorted_expected_df.to_csv(script_dir + "/test_output/sorted_expected.txt", sep="\t")
        
        tm.assert_frame_equal(sorted_expected_df, sorted_annot_df, check_names=False)
        
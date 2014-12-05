import ast
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
from jacquard.pivot_variants import PivotError, VariantPivoter, pivot, expand_format, merge_samples, create_initial_df, project_prepivot, build_pivoter, append_to_annot_df, melt_samples, validate_parameters, validate_format_tags, rearrange_columns, change_order, determine_input_keys, get_headers_and_readers

pd.set_option('chained_assignment', None)
TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))


def dataframe(input_data, sep="\t", index_col=None):
    def tupelizer(thing):
        if isinstance(thing, str) and thing.startswith("(") and thing.endswith(")"):
            return ast.literal_eval(thing)
        return thing

    df = pd.read_csv(StringIO(input_data), sep=sep, header=False, dtype='str', index_col=index_col)
    new_cols = [tupelizer(col) for col in list(df.columns.values)]
    df.columns = pd.core.index.Index(new_cols)
    
    return df
    
class VariantPivoterTestCase(unittest.TestCase):
    def test_build_transform(self):
        rows = ['COORDINATE']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        pivoter = VariantPivoter(rows, cols, pivot_values)
        
        input_string = \
'''COORDINATE	INFO	FORMAT	sample_A	sample_B
1	blah	DP:ESAF	10:0.2	100:0.2
2	blah	DP:ESAF	20:0.2	200:0.2
3	blah	DP:ESAF	30:0.2	300:0.2
4	blah	DP:ESAF	40:0.2	400:0.2'''
        df = dataframe(input_string)
        
        transform = pivoter._build_transform_method(rows, cols, pivot_values)
        actual_df = transform(df, "foo")
        actual_df.columns.names = [""]
        expected_string = \
'''COORDINATE	SAMPLE_NAME	DP
1	foo_sample_A	10
1	foo_sample_B	100
2	foo_sample_A	20
2	foo_sample_B	200
3	foo_sample_A	30
3	foo_sample_B	300
4	foo_sample_A	40
4	foo_sample_B	400'''
        expected_df = dataframe(expected_string)
        expected_df.columns.names = [""]

        tm.assert_frame_equal(expected_df, actual_df)

    def test_add_files(self):
        annot_df = pd.DataFrame()
        rows = ['COORDINATE']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        
        pivoter = VariantPivoter(rows, cols, pivot_values)
        sample_A_file = \
'''COORDINATE	FORMAT	Samp1
1	DP:ESAF	1:0.2
2	DP:ESAF	12:0.2
3	DP:ESAF	31:0.2
4	DP:ESAF	6:0.2'''
        sample_B_file = \
'''COORDINATE	FORMAT	Samp2
1	DP:ESAF	5:0.2
2	DP:ESAF	2:0.2
3	DP:ESAF	74:0.2
4	DP:ESAF	25:0.2'''

        pivoter.add_file("sample_A", StringIO(sample_A_file), 0)
        pivoter.add_file("sample_B", StringIO(sample_B_file), 0)
        
        actual_df = pivoter._combined_df
        actual_df.columns.names = [""]
        expected_string = \
'''COORDINATE	SAMPLE_NAME	DP
1	sample_A_Samp1	1
2	sample_A_Samp1	12
3	sample_A_Samp1	31
4	sample_A_Samp1	6
1	sample_B_Samp2	5
2	sample_B_Samp2	2
3	sample_B_Samp2	74
4	sample_B_Samp2	25
'''
        expected_df = dataframe(expected_string)
        expected_df.columns.names = [""]

        tm.assert_frame_equal(expected_df, actual_df)
    
    ##is_compatible
    def test_is_compatible_raiseIfMissingRequiredColumns(self): 
        rows = ['COORDINATE', 'foo']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        pivoter = VariantPivoter(rows, cols, pivot_values)
        
        input_string = \
'''COORDINATE	FORMAT	sample_A	sample_B
1	DP:ESAF	10:0.2	100:0.2
2	DP:ESAF	20:0.2	200:0.2
3	DP:ESAF	30:0.2	300:0.2
4	DP:ESAF	40:0.2	400:0.2'''
        df = dataframe(input_string)
        
        self.assertRaises(PivotError, pivoter.is_compatible, df, "foo")
        
    def test_check_required_columns_present(self):
        rows = ['Foo']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        pivoter = VariantPivoter(rows, cols, pivot_values)
       
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
        pivoter = VariantPivoter(rows, cols, pivot_values)
       
        expected_string = \
'''Foo	Bar	Baz	Blah
1	A	42	2
1	A	42	2'''

        df = dataframe(expected_string)

        self.assertRaises(PivotError, pivoter._check_pivot_is_unique, df)
    
    ##mult alts and mult genes
    def test_insert_mult_alts(self):
        rows = ['CHROM', 'POS', 'REF', 'ALT']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        
        input_string = \
'''CHROM	POS	REF	ALT	FORMAT	sample_A	sample_B
1	2	A	G	DP	134	57
1	3	C	G	DP	135	58
8	4	A	T	DP	136	59
13	5	T	AA	DP	137	60
13	5	T	AAAA	DP	137	60
'''
        df = dataframe(input_string)
        
        pivoter = VariantPivoter(rows, cols, pivot_values)
        actual_df = pivoter.insert_mult_alt_and_gene(df)
        
        expected_string = \
'''CHROM	POS	REF	ALT	Mult_Alt	FORMAT	sample_A	sample_B
1	2	A	G		DP	134	57
1	3	C	G		DP	135	58
8	4	A	T		DP	136	59
13	5	T	AA	True	DP	137	60
13	5	T	AAAA	True	DP	137	60
'''
        expected_df = dataframe(expected_string)
        expected_df.fillna(value="", inplace=True)

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_insert_mult_genes(self):
        rows = ['CHROM', 'POS', 'REF', 'ALT', 'GENE_SYMBOL']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        
        input_string = \
'''CHROM	POS	REF	ALT	GENE_SYMBOL	FORMAT	sample_A	sample_B
1	2	A	G	BRCA1	DP	134	57
1	3	C	G	EGFR	DP	135	58
8	4	A	T	TANK	DP	136	59
13	5	T	AA	ERBB2	DP	137	60
13	5	T	AA	CREBBP	DP	137	60
'''
        df = dataframe(input_string)
        
        pivoter = VariantPivoter(rows, cols, pivot_values)
        actual_df = pivoter.insert_mult_alt_and_gene(df)
        
        expected_string = \
'''CHROM	POS	REF	ALT	GENE_SYMBOL	Mult_Alt	Mult_Gene	FORMAT	sample_A	sample_B
1	2	A	G	BRCA1			DP	134	57
1	3	C	G	EGFR			DP	135	58
8	4	A	T	TANK			DP	136	59
13	5	T	AA	ERBB2		True	DP	137	60
13	5	T	AA	CREBBP		True	DP	137	60
'''
        expected_df = dataframe(expected_string)
        expected_df.fillna(value="", inplace=True)

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_insert_mult_alt_and_genes(self):
        rows = ['CHROM', 'POS', 'REF', 'ALT', 'GENE_SYMBOL']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        
        input_string = \
'''CHROM	POS	REF	ALT	GENE_SYMBOL	FORMAT	sample_A	sample_B
1	2	A	G	BRCA1	DP	134	57
1	3	C	G	EGFR	DP	135	58
8	4	A	T	TANK	DP	136	59
8	4	A	TTT	TANK	DP	136	59
13	5	T	AA	ERBB2	DP	137	60
13	5	T	AA	CREBBP	DP	137	60
'''
        df = dataframe(input_string)
        
        pivoter = VariantPivoter(rows, cols, pivot_values)
        actual_df = pivoter.insert_mult_alt_and_gene(df)
        
        expected_string = \
'''CHROM	POS	REF	ALT	GENE_SYMBOL	Mult_Alt	Mult_Gene	FORMAT	sample_A	sample_B
1	2	A	G	BRCA1			DP	134	57
1	3	C	G	EGFR			DP	135	58
8	4	A	T	TANK	True		DP	136	59
8	4	A	TTT	TANK	True		DP	136	59
13	5	T	AA	ERBB2		True	DP	137	60
13	5	T	AA	CREBBP		True	DP	137	60
'''
        expected_df = dataframe(expected_string)
        expected_df.fillna(value="", inplace=True)

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_label_mult_alt(self):
        rows = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
        cols = ["SAMPLE_NAME"]
        pivot_values = ['GT']
        pivoter = VariantPivoter(rows, cols, pivot_values)
    
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
1	2	A	T	foo	.	GT	0/1	10	100
1	2	A	G	foo	.	GT	0/1	10	100
2	3	A	T	foo	.	GT	0/1	20	200'''
        df = dataframe(input_string)
        
        cols_to_group = rows + ["WARNING/ERROR"]
        grouped = df.groupby(cols_to_group)
        df.set_index(cols_to_group, inplace=True)
        
        df = pivoter.label_mult(grouped, df, 'alt')
        df.fillna(value="", inplace=True)
        
        expected_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B	Mult_Alt
1	2	A	T	foo	.	GT	0/1	10	100	True
1	2	A	G	foo	.	GT	0/1	10	100	True
2	3	A	T	foo	.	GT	0/1	20	200	'''
        expected_df = dataframe(expected_string)
        expected_df.set_index(cols_to_group, inplace=True)
        expected_df.fillna(value="", inplace=True)
        
        tm.assert_frame_equal(expected_df, df)
    
    def test_label_mult_gene(self):
        rows = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
        cols = ["SAMPLE_NAME"]
        pivot_values = ['GT']
        pivoter = VariantPivoter(rows, cols, pivot_values)
    
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
1	2	A	T	foo	.	GT	0/1	10	100
1	2	A	T	bar	.	GT	0/1	10	100
2	3	A	T	foo	.	GT	0/1	20	200'''
        df = dataframe(input_string)
        
        cols_to_group = rows + ["WARNING/ERROR"]
        grouped = df.groupby(cols_to_group)
        df.set_index(cols_to_group, inplace=True)
        
        df = pivoter.label_mult(grouped, df, 'gene')
        df.fillna(value="", inplace=True)
        
        expected_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B	Mult_Gene
1	2	A	T	foo	.	GT	0/1	10	100	True
1	2	A	T	bar	.	GT	0/1	10	100	True
2	3	A	T	foo	.	GT	0/1	20	200	'''
        expected_df = dataframe(expected_string)
        expected_df.set_index(cols_to_group, inplace=True)
        expected_df.fillna(value="", inplace=True)
        
        tm.assert_frame_equal(expected_df, df)
    
    def test_create_mult_dict_alt(self):
        rows = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
        cols = ["SAMPLE_NAME"]
        pivot_values = ['GT']
        pivoter = VariantPivoter(rows, cols, pivot_values)
    
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
sample1	1	2	A	T	foo	.	GT	0/1	10	100
sample1	1	2	A	G	foo	.	GT	0/1	10	100
sample2	2	3	A	T	foo	.	GT	0/1	20	200'''
        df = dataframe(input_string)
        
        cols_to_group = rows + ["WARNING/ERROR"]
        grouped = df.groupby(cols_to_group)
        
        val_index = 3
        preliminary_dict = {}
        mult_dict = {}
        actual_dict = pivoter._create_mult_dict(grouped, val_index, preliminary_dict, mult_dict)
        
        expected_dict = {"1|2|A|foo|." : ["G", "T"]}
        
        self.assertEquals(expected_dict, actual_dict)
        
    def test_create_mult_dict_gene(self):
        rows = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
        cols = ["SAMPLE_NAME"]
        pivot_values = ['GT']
        pivoter = VariantPivoter(rows, cols, pivot_values)
    
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
sample1	1	2	A	T	foo	.	GT	0/1	10	100
sample1	1	2	A	T	bar	.	GT	0/1	10	100
sample2	2	3	A	T	foo	.	GT	0/1	20	200'''
        df = dataframe(input_string)
        
        cols_to_group = rows + ["WARNING/ERROR"]
        grouped = df.groupby(cols_to_group)
        
        val_index = 4
        preliminary_dict = {}
        mult_dict = {}
        actual_dict = pivoter._create_mult_dict(grouped, val_index, preliminary_dict, mult_dict)
        
        expected_dict = {"1|2|A|T|." : ["bar", "foo"]}
        
        self.assertEquals(expected_dict, actual_dict)
    
    def test_determine_row_key(self):
        rows = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
        cols = ["SAMPLE_NAME"]
        pivot_values = ['GT']
        pivoter = VariantPivoter(rows, cols, pivot_values)
        
        row_key = ["1", "2", "A", "T", "."]
        val = "bar"
        val_index = 4
        actual_key = pivoter._determine_row_key(row_key, val, val_index)
        
        expected_key = ("1", "2", "A", "T", "bar", ".")
        
        self.assertEquals(expected_key, actual_key)
        
    ##validate data
    def test_validate_annotations(self):
        rows = ['COORDINATE']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        
        input_string = \
'''COORDINATE	INFO	FORMAT	sample_A	sample_B
1	blah	DP	10	100
1	blah	DP	10	100
2	blah	DP	20	200
3	blah	DP	30	300
4	blah	DP	40	400'''
        df = dataframe(input_string)
        
        pivoter = VariantPivoter(rows, cols, pivot_values, combined_df=pd.DataFrame(), annot_df=df)
        actual_df = pivoter.validate_annotations()
        
        expected_string = \
'''COORDINATE	INFO	FORMAT	sample_A	sample_B
1	blah	DP	10	100
2	blah	DP	20	200
3	blah	DP	30	300
4	blah	DP	40	400'''
        expected_df = dataframe(expected_string)
        
        tm.assert_frame_equal(expected_df, actual_df)
    
    #memory error for large files
    def test_validate_sample_data_non_unique_rows(self):
        rows = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
        cols = ["SAMPLE_NAME"]
        pivot_values = ['GT']
        
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	SAMPLE_NAME	SAMPLE_DATA
1	2	A	CG	foo	.	GT	sampleA	10
1	2	A	T	foo	.	GT	sampleA	13
1	2	A	T	foo	.	GT	sampleA	12
2	3	C	G	foo	.	GT	sampleB	11'''
        df = dataframe(input_string)
        combined_df = df
        
        pivoter = VariantPivoter(rows, cols, pivot_values, combined_df)
        actual_df = pivoter.validate_sample_data()
        
        expected_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	SAMPLE_NAME	SAMPLE_DATA
1	2	A	CG	foo	.	GT	sampleA	10
1	2	A	T	foo	.	GT	sampleA	^
2	3	C	G	foo	.	GT	sampleB	11'''
        expected_df = dataframe(expected_string)
        print actual_df
        print expected_df
        tm.assert_frame_equal(expected_df, actual_df)
        
    # not sure how to test if numpy.ndarray
    def Ytest_validate_sample_data_non_unique_cells(self):
        rows = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"]
        cols = ["SAMPLE_NAME"]
        pivot_values = ['GT']
          
        input_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	SAMPLE_NAME	SAMPLE_DATA
1	2	A	T	foo	.	GT	sampleA	{0}
1	2	A	C	foo	.	GT	sampleA	11
2	3	C	G	foo	.	GT	sampleB	10'''.format(np.ndarray(shape=(1,2)))
        df = dataframe(input_string)
          
        pivoter = VariantPivoter(rows, cols, pivot_values, df)
        actual_df = pivoter.validate_sample_data()
          
        expected_string = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	SAMPLE_NAME	SAMPLE_DATA
1	2	A	T	foo	.	GT	sampleA	^
1	2	A	C	foo	.	GT	sampleA	11
2	3	C	G	foo	.	GT	sampleB	11'''
        expected_df = dataframe(expected_string)
        print actual_df.applymap(lambda x: type(x))
          
        print actual_df
        tm.assert_frame_equal(expected_df, actual_df)
    
class PivotTestCase(unittest.TestCase):
    ##validate parameters
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

    def test_validate_format_tags_invalid(self):
        first_line = ["1\t24\tA\tC\tGT:DP\tfoo;bar\t1/1:258"]
        pivot_values = ["foo"]
        fields = ["CHROM", "POS", "REF", "ALT", "FORMAT", "INFO", "sample2"]
        actual_invalid_tags = validate_format_tags(first_line, pivot_values, fields)
        
        expected_invalid_tags = ["foo"]
        
        self.assertEqual(expected_invalid_tags, actual_invalid_tags)
        
    def test_validate_format_tags_valid(self):
        first_line = ["1\t24\tA\tC\tGT:DP\tfoo;bar\t1/1:258"]
        pivot_values = ["GT"]
        fields = ["CHROM", "POS", "REF", "ALT", "FORMAT", "INFO", "sample2"]
        actual_invalid_tags = validate_format_tags(first_line, pivot_values, fields)
        
        expected_invalid_tags = []
        
        self.assertEqual(expected_invalid_tags, actual_invalid_tags)
    
    ##determine input keys
    def test_determine_input_keys_txt(self):
        input_dir = TEST_DIRECTORY + "/reference_files/test_input/test_input_keys_txt"
        actual_lst = determine_input_keys(input_dir)
        
        expected_lst = ["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL", "SnpEff_WARNING/ERROR"]
        
        self.assertEquals(expected_lst, actual_lst)
        
    def test_determine_input_keys_vcf(self):
        input_dir = TEST_DIRECTORY + "/reference_files/test_input/test_input_keys_vcf"
        actual_lst = determine_input_keys(input_dir)
        
        expected_lst = ["CHROM", "POS", "REF", "ALT"]
        
        self.assertEquals(expected_lst, actual_lst)
        
    def test_determine_input_keys_invalid(self):
        input_dir = TEST_DIRECTORY + "/reference_files/test_input/test_input_keys_invalid"
        
        self.assertRaises(PivotError, determine_input_keys, input_dir)
    
    ##get headers, readers
    def test_get_headers_and_readers(self):
        input_dir = TEST_DIRECTORY + "/reference_files/test_input/test_input_keys_txt"
        sample_file_readers, headers, header_names, first_line = get_headers_and_readers(input_dir)
        
        self.assertEquals([input_dir + "/foo1.txt", input_dir + "/foo2.txt"], sample_file_readers)
        self.assertEquals([3,2], headers)
        self.assertEquals("CHROM	POS	REF	ALT	GENE_SYMBOL	FORMAT	Sample_2384	Sample_2385", header_names.rstrip())
        self.assertEquals("1	2342	A	T	EGFR	GT:DP	1/1:241	0/1:70", first_line[0].rstrip())
        self.assertEquals("1	134	G	C	EGFR	GT:DP	1/1:242	0/1:546", first_line[1].rstrip())
    
    def test_build_pivoter_invalidHeaderRaisesPivotError(self):
        input_string = \
'''COORDINATE	FORMAT	sample_A	sample_B
1	GT:ESAF	10:0.2	100:0.2
2	GT:ESAF	20:0.2	200:0.2
3	GT:ESAF	30:0.2	300:0.2
4	GT:ESAF	40:0.2	400:0.2'''
        
        input_df = dataframe(input_string)
        input_keys = ['CHROM', 'POS']
        pivot_values = ["GT"]
            
        self.assertRaises(PivotError, build_pivoter, input_keys, pivot_values, input_df)

    ##create_initial_df
    def test_build_pivoter_keepWarnings(self):
        input_string = \
'''COORDINATE	FORMAT	sample_A	sample_B	WARNING/ERROR
1	GT:ESAF	10:0.2	100:0.2	.
2	GT:ESAF	20:0.2	200:0.2	warning error line
3	GT:ESAF	30:0.2	300:0.2	.
4	GT:ESAF	40:0.2	400:0.2	.'''
        input_df = dataframe(input_string)
        expected_string = \
'''COORDINATE	FORMAT	_sample_A	_sample_B	SnpEff_WARNING/ERROR
1	GT:ESAF	10:0.2	100:0.2	.
2	GT:ESAF	20:0.2	200:0.2	warning error line
3	GT:ESAF	30:0.2	300:0.2	.
4	GT:ESAF	40:0.2	400:0.2	.'''
        expected_df = dataframe(expected_string)
        input_keys = ["COORDINATE"]
        pivot_values = ["GT"]
        pivoter,actual_df = build_pivoter(input_keys, pivot_values, input_df)

        tm.assert_frame_equal(expected_df, actual_df)
        
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
1	2	GT:DP:ESAF	A:1:0.2
1	3	GT:DP:ESAF	B:2:0.2
1	5	GT:DP:ESAF	:3:0.2
13	3	GT:DP:ESAF	D::0.2'''
        input_df = pd.read_csv(StringIO(dataString), sep="\t", header=False)
        rows=["CHROM", "POS"]
        actual_df = expand_format(input_df, ["GT"], rows, "foo")
 
        expected_string = \
'''CHROM	POS	SAMPLE_NAME	DP	ESAF	GT
1	2	foo_Sample1	1	0.2	A
1	3	foo_Sample1	2	0.2	B
1	5	foo_Sample1	3	0.2	
13	3	foo_Sample1		0.2	D'''

        expected_df = pd.read_csv(StringIO(expected_string), sep="\t", header=False, dtype={'DP':str, 'ESAF':str})
        expected_df.fillna(value="", inplace=True)

        tm.assert_series_equal(expected_df.ix[:,0], actual_df.ix[:,0])
        tm.assert_series_equal(expected_df.ix[:,1], actual_df.ix[:,1])
        tm.assert_series_equal(expected_df.ix[:,2], actual_df.ix[:,2])
        tm.assert_series_equal(expected_df.ix[:,3], actual_df.ix[:,3])
        tm.assert_series_equal(expected_df.ix[:,4], actual_df.ix[:,4])
        tm.assert_series_equal(expected_df.ix[:,5], actual_df.ix[:,5])
        
    def test_expand_format_multipleFormats(self):
        dataString = \
'''CHROM	POS	FORMAT	Sample1
1	1	GT:DP:ESAF	A:1:0.2
1	2	GT:DP:ESAF	B:2:0.2
1	3	GT:DP:ESAF	:3:0.2
1	4	GT:DP:ESAF	D::0.2'''
        input_df = pd.read_csv(StringIO(dataString), sep="\t", header=False)
        rows = ["CHROM", "POS"]
        actual_df = expand_format(input_df, ["GT", "DP"], rows, "foo")

        expectedString = \
'''CHROM\tPOS\tSAMPLE_NAME\tDP\tESAF\tGT
1	1	foo_Sample1	1	0.2	A
1	2	foo_Sample1	2	0.2	B
1	3	foo_Sample1	3	0.2	
1	4	foo_Sample1		0.2	D'''
        expected_df = pd.read_csv(StringIO(expectedString), sep="\t", header=False, dtype={'DP':str, 'ESAF':str})
        expected_df.fillna(value="", inplace=True)

        tm.assert_series_equal(expected_df.ix[:,0], actual_df.ix[:,0])
        tm.assert_series_equal(expected_df.ix[:,1], actual_df.ix[:,1])
        tm.assert_series_equal(expected_df.ix[:,2], actual_df.ix[:,2])
        tm.assert_series_equal(expected_df.ix[:,3], actual_df.ix[:,3])
        tm.assert_series_equal(expected_df.ix[:,4], actual_df.ix[:,4])
        tm.assert_series_equal(expected_df.ix[:,5], actual_df.ix[:,5])
        
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
chr1	1	A	T	GT:ESAF	0/1:0.2	1/1:0.2
chr1	2	A	T	GT:ESAF	0/1:0.2	0/1:0.2
chr1	3	A	T	GT:ESAF	1/1:0.2	1/1:0.2
chr1	4	A	T	GT:ESAF	1/1:0.2	0/1:0.2'''
        df = pd.read_csv(StringIO(dataString), sep="\t", header=False, dtype='str')
        
        actual_df = melt_samples(df, "foo")
        
        expected_dataString = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	FORMAT	SAMPLE_NAME	SAMPLE_DATA
chr1	1	A	T	GT:ESAF	foo_Sample_1	0/1:0.2
chr1	2	A	T	GT:ESAF	foo_Sample_1	0/1:0.2
chr1	3	A	T	GT:ESAF	foo_Sample_1	1/1:0.2
chr1	4	A	T	GT:ESAF	foo_Sample_1	1/1:0.2
chr1	1	A	T	GT:ESAF	foo_Sample_2	1/1:0.2
chr1	2	A	T	GT:ESAF	foo_Sample_2	0/1:0.2
chr1	3	A	T	GT:ESAF	foo_Sample_2	1/1:0.2
chr1	4	A	T	GT:ESAF	foo_Sample_2	0/1:0.2'''
        expected_df = pd.read_csv(StringIO(expected_dataString), sep="\t", header=False, dtype='str')

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_melt_samples_trailing_field(self):
        dataString = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	FORMAT	Sample_1	Sample_2	Sample
chr1	1	A	T	GT:DP:ESAF	0/1:2:0.2	1/1:12:0.2	foo
chr1	2	A	T	GT:DP:ESAF	0/1:3:0.2	0/1:13:0.2	foo
chr1	3	A	T	GT:DP:ESAF	1/1:4:0.2	1/1:14:0.2	foo
chr1	4	A	T	GT:DP:ESAF	1/1:5:0.2	0/1:15:0.2	foo'''
        df = pd.read_csv(StringIO(dataString), sep="\t", header=False, dtype='str')
        
        actual_df = melt_samples(df, "foo")
        
        expected_dataString = \
'''CHROM	POS	REF	ANNOTATED_ALLELE	FORMAT	Sample	SAMPLE_NAME	SAMPLE_DATA
chr1	1	A	T	GT:DP:ESAF	foo	foo_Sample_1	0/1:2:0.2
chr1	2	A	T	GT:DP:ESAF	foo	foo_Sample_1	0/1:3:0.2
chr1	3	A	T	GT:DP:ESAF	foo	foo_Sample_1	1/1:4:0.2
chr1	4	A	T	GT:DP:ESAF	foo	foo_Sample_1	1/1:5:0.2
chr1	1	A	T	GT:DP:ESAF	foo	foo_Sample_2	1/1:12:0.2
chr1	2	A	T	GT:DP:ESAF	foo	foo_Sample_2	0/1:13:0.2
chr1	3	A	T	GT:DP:ESAF	foo	foo_Sample_2	1/1:14:0.2
chr1	4	A	T	GT:DP:ESAF	foo	foo_Sample_2	0/1:15:0.2'''
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
    
    ##rearrange columns
    def test_rearrange_columns(self):
        input_string = \
'''CHROM	POS	REF	ALT	FORMAT	('DP', 'sample_A')	('DP', 'sample_B')
1	2	A	G	DP	134	57
1	3	C	G	DP	135	58
8	4	A	T	DP	136	59
13	5	T	AA	DP	137	60
13	5	T	AAAA	DP	137	60
'''
        df = dataframe(input_string)
        
        actual_df = rearrange_columns(df)
        
        expected_string = \
'''CHROM	POS	REF	ALT	DP_sample_A	DP_sample_B	FORMAT
1	2	A	G	134	57	DP
1	3	C	G	135	58	DP
8	4	A	T	136	59	DP
13	5	T	AA	137	60	DP
13	5	T	AAAA	137	60	DP
'''
        expected_df = dataframe(expected_string)

        tm.assert_frame_equal(expected_df, actual_df)
    
    def test_change_order(self):
        lst = ["CHROM", "POS", "REF", "ALT", "FORMAT", "INFO", "DP_fname1_sample1", "DP_fname2_sample2"]
        pivot_columns = ["DP_fname1_sample1", "DP_fname2_sample2"]
        actual_lst = change_order(lst, pivot_columns)
        
        expected_lst = ["CHROM", "POS", "REF", "ALT", "DP_fname1_sample1", "DP_fname2_sample2", "FORMAT", "INFO"]
        
        self.assertEqual(expected_lst, actual_lst)
        
    def test_change_order_formatsInFilenamesOkay(self):
        lst = ["CHROM", "POS", "REF", "ALT", "FORMAT", "INFO", "DP_fnameDPGQ_sample1", "GQ_fnameDPGQ_sample1"]
        pivot_columns = ["DP_fnameDPGQ_sample1", "GQ_fnameDPGQ_sample1"]
        actual_lst = change_order(lst, pivot_columns)
        
        expected_lst = ["CHROM", "POS", "REF", "ALT", "DP_fnameDPGQ_sample1", "GQ_fnameDPGQ_sample1", "FORMAT", "INFO"]
        
        self.assertEqual(expected_lst, actual_lst)
        
    ##annot df
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
        
        file_list = [script_dir + "/reference_files/test_input/P2_test_input.txt", script_dir + "/reference_files/test_input/P5_test_input.txt"]
        for file in file_list:
            df = pd.read_csv(file, sep="\t", header=1, dtype='str', index_col=False)
            annot_df = append_to_annot_df(df, annot_df)

        expected_df = pd.read_csv(script_dir + "/reference_files/test_annotation/P2P5_combined_annotation.txt", sep="\t", header=False, dtype='str', index_col=False)

        try:
            expected_df["CHROM"] = expected_df["CHROM"].apply(lambda x: x.replace("chr", ""))
            annot_df["CHROM"] = annot_df["CHROM"].apply(lambda x: x.replace("chr", ""))
           
            expected_df["CHROM"] = expected_df["CHROM"].apply(lambda x: int(x))
            annot_df["CHROM"] = annot_df["CHROM"].apply(lambda x: int(x))
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

        sorted_annot_df.to_csv(script_dir + "/reference_files/test_output/sorted_annot.txt", sep="\t")
        sorted_expected_df.to_csv(script_dir + "/reference_files/test_output/sorted_expected.txt", sep="\t")
        
        tm.assert_frame_equal(sorted_expected_df, sorted_annot_df, check_names=False)
        
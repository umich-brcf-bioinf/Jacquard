#!/usr/bin/env python
import pandas as pd
from pandas import *
import unittest
from pandas.util.testing import assert_frame_equal
import pandas.util.testing as tm
from bin.pivot import PivotError, VariantPivoter, EpeeVariantPivoter, VcfVariantPivoter, pivot, expand_format, merge_samples, create_initial_df, project_prepivot, build_pivoter
from StringIO import StringIO
import pprint


def dataframe(input_data, sep="\t", index_col=None):
    return pd.read_csv(StringIO(input_data), sep=sep, header=False, dtype='str', index_col=index_col)

class VariantPivoterTestCase(unittest.TestCase):
    def test_add_files(self):
        rows = ['COORDINATE']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        transform = lambda x : x
        pivoter = VariantPivoter(rows, cols, pivot_values, transform)
        sample_A_file = \
'''COORDINATE	DP
1	10
2	20
3	30
4	40'''
        sample_B_file = \
'''COORDINATE	DP
1	100
2	200
3	300
4	400'''

        pivoter.add_file("sample_A", StringIO(sample_A_file))
        pivoter.add_file("sample_B", StringIO(sample_B_file))
        actual_df = pivoter.pivot()

        actual_df = actual_df['DP']
        actual_stringIO = StringIO()
        actual_df.to_csv(actual_stringIO, sep="\t")
        actual_flattened_df = pd.read_csv(StringIO(actual_stringIO.getvalue()), header=False, sep="\t", dtype='str')

        expected_string = \
'''COORDINATE	sample_A	sample_B
1	10	100
2	20	200
3	30	300
4	40	400'''
        tm.assert_frame_equal(dataframe(expected_string), actual_flattened_df)
        
    def test_add_files_transformAsFilter(self):
        rows = ['COORDINATE']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        transform = lambda df : df[df["FILTER_COL"]=='.']
        pivoter = VariantPivoter(rows, cols, pivot_values, transform)
        sample_A_file = \
'''COORDINATE	DP	FILTER_COL
1	10	.
2	20	.
3	30	Yep
4	40	.
'''
        sample_B_file = \
'''COORDINATE	DP	FILTER_COL
1	100	.
2	200	Yep
3	300	.
4	400	.
'''
        pivoter.add_file("sample_A", StringIO(sample_A_file))
        pivoter.add_file("sample_B", StringIO(sample_B_file))
        actual_df = pivoter.pivot()

        actual_df = actual_df['DP']
        actual_stringIO = StringIO()
        actual_df.to_csv(actual_stringIO, sep="\t")
        actual_flattened_df = pd.read_csv(StringIO(actual_stringIO.getvalue()), header=False, sep="\t", dtype='str')

        expected_string = \
'''COORDINATE	sample_A	sample_B
1	10	100
2	20	NaN
3	NaN	300
4	40	400'''
        tm.assert_frame_equal(dataframe(expected_string), actual_flattened_df)

       
    def test_is_compatible_raiseIfMissingRequiredColumns(self): 
        rows = ['Foo']
        cols = ['SAMPLE_NAME']
        pivot_values = ['DP']
        transform = lambda x : x
        pivoter = VariantPivoter(rows, cols, pivot_values, transform)
        
        input_string = \
'''COORDINATE	sample_A	sample_B
1	10	100
2	20	200
3	30	300
4	40	400'''
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
        

class EpeeVariantPivoterTestCase(unittest.TestCase):

    def test_build_transform_method(self):
        input_df = dataframe(
            "SAMPLE_NAME|#CHROM|POS|WARNING/ERROR|FORMAT|Sample1\n" + \
            "SampleA|1|41|warn|GT:DP|0/41:410\n" + \
            "SampleB|1|42|.|GT:DP|0/42:420\n",'|')
        transform = EpeeVariantPivoter._build_transform_method(['#CHROM','POS'], ['SAMPLE_NAME'],['DP'])

        actual_df = transform(input_df)
        
        expected_df = dataframe(
            "|SAMPLE_NAME|#CHROM|POS|DP\n" + \
            "1|SampleB|1|42|420",sep='|', index_col=0)
        tm.assert_frame_equal(expected_df, actual_df)


    def test_exclude_errors_warnings_throwsIfColumnMissing(self):
        input_df = dataframe(
            "#CHROM|POS\n" + \
            "1|A\n" + \
            "2|B",'|')
        filter = EpeeVariantPivoter(["GT"])._pivoter._transform

        self.assertRaises(PivotError, filter, input_df)

    def test_is_compatible_trueIfCompatible(self):
        pivot_values = ['DP']
        pivoter = EpeeVariantPivoter(pivot_values)        
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
sample1	1	2	A	T	foo	.	GT	0/1	10	100
sample2	1	2	A	T	foo	WARN	GT	0/1	10	100
sample3	2	3	A	T	foo	.	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(True, pivoter.is_compatible(df))
    
    def test_is_compatible_falseIfMissingColumns(self):
        pivot_values = ['DP']
        pivoter = EpeeVariantPivoter(pivot_values)        
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
sample1	1	2	A	foo	.	GT	0/1	10	100
sample2	1	2	A	foo	WARN	GT	0/1	10	100
sample3	2	3	A	foo	.	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(False, pivoter.is_compatible(df))
        
    def test_is_compatible_trueIfWarningDuplicates(self):
        pivot_values = ['DP']
        pivoter = EpeeVariantPivoter(pivot_values)        
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
sample1	1	2	A	T	foo	.	GT	0/1	10	100
sample1	1	2	A	T	foo	WARN	GT	0/1	10	100
sample3	2	3	A	T	foo	.	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(True, pivoter.is_compatible(df))
        
    def test_is_compatible_falseIfDuplicates(self):
        pivot_values = ['DP']
        pivoter = EpeeVariantPivoter(pivot_values)        
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
sample1	1	2	A	T	foo	.	GT	0/1	10	100
sample1	1	2	A	T	foo	.	GT	0/1	10	100
sample3	2	3	A	T	foo	.	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(False, pivoter.is_compatible(df))

        
class VcfVariantPivoterTestCase(unittest.TestCase):
    def test_is_compatible_trueIfCompatible(self):
        pivot_values = ['DP']
        pivoter = VcfVariantPivoter(pivot_values)        
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ALT	GENE_SYMBOL	FORMAT	Sample1	sample_A	sample_B
sample1	1	2	A	T	foo	GT	0/1	10	100
sample2	1	2	A	T	foo	GT	0/1	10	100
sample3	2	3	A	T	foo	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(True, pivoter.is_compatible(df))  
        
    def test_is_compatible_falseIfDuplicates(self):
        pivot_values = ['DP']
        pivoter = EpeeVariantPivoter(pivot_values)        
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	FORMAT	Sample1	sample_A	sample_B
sample1	1	2	A	T	foo	GT	0/1	10	100
sample1	1	2	A	T	foo	GT	0/1	10	100
sample3	2	3	A	T	foo	GT	0/1	20	200
'''
        df = dataframe(input_string)

        self.assertEqual(False, pivoter.is_compatible(df))
        
class PivotTestCase(unittest.TestCase):
    def test_build_pivoter_invalidHeaderRaisesPivotError(self):
        input_string = \
'''COORDINATE	sample_A	sample_B
1	10	100
2	20	200
3	30	300
4	40	400'''
            
        pivot_values = ["GT"]
            
        self.assertRaises(PivotError, build_pivoter, "sampleA", StringIO(input_string), pivot_values)

    def test_build_pivoter_epeeHeader(self):
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	WARNING/ERROR	FORMAT	Sample1	sample_A	sample_B
samp1	1	2	A	T	foo	WARN	GT	0/1	10	100
samp2	2	3	A	T	foo	.	GT	0/1	20	200
samp3	3	4	A	T	foo	.	GT	0/1	30	300
samp4	4	5	A	T	foo	.	GT	0/1	40	400'''
       
        pivot_values = ["GT"]
        
        pivoter = build_pivoter("sampleA", StringIO(input_string), pivot_values)

        self.assertIsInstance(pivoter, EpeeVariantPivoter)
        
    def test_build_pivoter_vcfHeader(self):
        input_string = \
'''CHROM	POS	REF	ALT	sample_A	sample_B
1	2	A	T	10	100
2	3	A	T	20	200
3	4	A	T	30	300
4	5	A	T	40	400'''
        pivot_values = ["GT"]
        
        pivoter = build_pivoter("sampleA", StringIO(input_string), pivot_values)
        
        self.assertIsInstance(pivoter, VcfVariantPivoter)
            
        
    ##create_initial_df
    def test_create_initial_df(self):
        reader = StringIO( 
'''#CHROM	POS	REF
1	42	A
2	43	C''');
        
        path = 'foo/bar/sampleA_AllGenes.txt'
        actual_df = create_initial_df(path, reader)
        
        expected_data = StringIO(
'''#CHROM	POS	REF	SAMPLE_NAME
1	42	A	sampleA_AllGenes.txt
2	43	C	sampleA_AllGenes.txt''');

        expected_df = pd.read_csv(expected_data, sep="\t", header=False, dtype='str')
        print expected_df
        print actual_df
        tm.assert_frame_equal(expected_df, actual_df)

    ##expand_format   
    def test_expand_format(self):
        dataString = \
'''SAMPLE_NAME	FORMAT	Sample1
sample1	GT:DP	A:1
sample2	GT:DP	B:2
sample3	DP	3
sample6	GT	D'''
        input_df = pd.read_csv(StringIO(dataString), sep="\t", header=False)

        actual_df = expand_format(input_df, "FORMAT", "Sample1", ["GT"])
 
        expected_string = \
'''\tSAMPLE_NAME\tFORMAT\tSample1\tGT
0	sample1	GT:DP	A:1	A
1	sample2	GT:DP	B:2	B
2	sample3	DP	3	NaN
3	sample6	GT	D	D'''
        expected_df = dataframe(expected_string, index_col=0)        
        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_expand_format_multipleFormats(self):
        dataString = \
'''SAMPLE_NAME	FORMAT	Sample1
sample1	GT:DP	A:1
sample2	GT:DP	B:2
sample3	DP	3
sample6	GT	D'''
        input_df = pd.read_csv(StringIO(dataString), sep="\t", header=False)

        actual_df = expand_format(input_df, "FORMAT", "Sample1", ["GT", "DP"])

        expectedString = \
'''\tSAMPLE_NAME\tFORMAT\tSample1\tGT\tDP
0	sample1	GT:DP	A:1	A	1
1	sample2	GT:DP	B:2	B	2
2	sample3	DP	3	NaN	3
3	sample6	GT	D	D	NaN'''
        expected_df = pd.read_csv(StringIO(expectedString), sep="\t", header=False, index_col=0, dtype={'DP':str})        
        tm.assert_frame_equal(expected_df, actual_df)

    def test_expand_format_invalidFormatsIgnored(self):
        dataString = \
'''SAMPLE_NAME	FORMAT	Sample1
sample1	GT:DP	A:1'''
        input_df = pd.read_csv(StringIO(dataString), sep="\t", header=False)

        actual_df = expand_format(input_df, "FORMAT", "Sample1", ["GT", "DP", "FOO"])

        expectedString = \
'''\tSAMPLE_NAME\tFORMAT\tSample1\tGT\tDP
0	sample1	GT:DP	A:1	A	1'''
        expected_df = pd.read_csv(StringIO(expectedString), sep="\t", header=False, index_col=0, dtype={'DP':str})        
        tm.assert_frame_equal(expected_df, actual_df)

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
        reduced_df = project_prepivot(df, pivot_values, rows, columns)
        
        expected_dataString = \
'''SAMPLE_NAME	CHROM	POS	REF	ALT	FORMAT	GT
sample1	chr1	1	A	T	GT:DP	0/1
sample2	chr1	2	A	T	GT:DP	0/1
sample3	chr1	3	A	T	DP	1/1
sample6	chr1	4	A	T	GT	1/1'''
        expected_df = pd.read_csv(StringIO(expected_dataString), sep="\t", header=False)

        
        
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

    def test_merge_samples_PivotErrorIfColumnsDoNotMatch(self): #########
        file_readers = {
            'foo/sampleA':StringIO("#CHROM	COL_A\n1	A"),
            'foo/sampleB':StringIO("#CHROM	COL_B\n1	B")
            }

        # self.assertRaises(PivotError, merge_samples, file_readers)

        
    def test_merge_samples_raisesOnDifferentNumberOfColumns(self): #########
        sampleA_data = StringIO( 
'''#CHROM	POS	REF	FORMAT	SAMPLE_A
1	42	A	GT	0/1
2	43	C	GT	1/1''');
        sampleB_data = StringIO( 
'''#CHROM	POS	REF	FORMAT	SAMPLE_A	SAMPLE_B
1	44	G	GT	1/1	1/1
Y	45	T	GT	0/1	1/1''');
        
        file_readers = {
            'foo/bar/sampleA_AllGenes.txt':sampleA_data,
            'foo/bar/sampleB_AllGenes.txt':sampleB_data
            }

        # self.assertRaises(PivotError, merge_samples, file_readers)
        # self.assertEqual(1,2)
   
        
    def Xtest_pivot(self):
        dataString = \
'''SAMPLE_NAME	CHROM	POS	REF	ALT	FORMAT	Sample1	GT	DP
sample1	chr1	2	A	C	GT:DP	A:1	A	1
sample1	chr2	3	T	C	GT:DP	B:2	B	2
sample2	chr1	2	A	C	GT:DP	C:3	C	3
sample2	chr2	3	T	C	GT:DP	D:4	D	4'''
        input_df = pd.read_csv(StringIO(dataString), sep="\t", header=False)

        actual_df = pivot(input_df, 0, ["GT"])
        #actual_df = actual_df["GT"]
        print actual_df
        pp = pprint.PrettyPrinter(depth=6)
        pp.pprint(actual_df.to_dict())
        expectedString = \
'''SAMPLE_NAME				sample1	sample2
CHROM	POS	REF	ALT		
chr1	2	A	C	A	C
chr2	3	T	C	B	D'''
        expected_df = pd.read_csv(StringIO(expectedString), sep="\t", header=[0,1], index_col=0, dtype={'DP':str})
        expected_dict = {
            ('GT', 'sample1'): {('chr1', 2.0, 'A', 'C'): 'A',('chr2', 3.0, 'T', 'C'): 'B'},
            ('GT', 'sample2'): {('chr1', 2.0, 'A', 'C'): 'C',('chr2', 3.0, 'T', 'C'): 'D'}}
        expected_df = DataFrame.from_dict(expected_dict)
        print expected_df
        self.assertEqual(1,2)

       
        tm.assert_frame_equal(expected_df, actual_df)

        
        
    def Xtest_pivot(self):
        data = {
            'SAMPLE_NAME' : ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6'],
            'CHROM' : ['chr1', 'chr1', 'chr2', 'chr14', 'chr5', 'chr38'],
            'POS' : ['1', '2', '3', '3', '4', '5'],
            'REF' : ['A', 'T', 'G', 'C', 'A', 'G'],
            'ALT' : ['T', 'GC', 'C', 'A', 'T', 'C'],
            'WARNING/ERROR' : ['A', '.', 'C', '.', '.', '.'],
            'FORMAT' : ['GT:DP', 'GT:DP', 'GT:DP', 'GT:DP', 'GT:DP', 'GT:DP'],
            'Sample1' : ['0/1:23', '0/1:24', '1/1:25', '0/1:26', '1/1:27', '0/1:28'],
        }
        df = DataFrame(data) 
        # print df
        # exit(1)
        epee = 0
        format = ["GT"]
        actual_pivoted = pivot(df, epee, format)
        
        # for row, col in df.T.iteritems():
            # df.loc[row, "GT"] = 
        # expected_pivoted = df.pivot_table(df, values="GT", rows=["CHROM", "POS", "REF", "ANNOTATED_ALLELE", "GENE_SYMBOL"], cols=["SAMPLE_NAME"], aggfunc=lambda x: x)
        
        expected_data = {
            'SAMPLE_NAME' : ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6'],
            'CHROM' : ['chr1', 'chr1', 'chr2', 'chr14', 'chr5', 'chr38'],
            'POS' : ['1', '2', '3', '3', '4', '5'],
            'REF' : ['A', 'T', 'G', 'C', 'A', 'G'],
            'GT' : ['0/1', '0/1', '1/1', '0/1', '1/1', '0/1'],
        }
        expected_df = DataFrame(expected_data)
        # expected_pivoted = expected_df.set_index(["CHROM", "POS", "REF"])
        # expected_pivoted = expected_df.pivot(index=["CHROM", "POS", "REF"], columns=["SAMPLE_NAME"], values=["GT"])
                
        # print actual_pivoted
        # print expected_df
        # tm.assert_frame_equal(expected_pivoted, actual_pivoted)
    

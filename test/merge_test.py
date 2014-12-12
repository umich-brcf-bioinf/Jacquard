from argparse import Namespace
from collections import OrderedDict
import glob
import os
import pandas as pd
import pandas.util.testing as tm
from StringIO import StringIO
import sys
from testfixtures import TempDirectory
import unittest

from jacquard.merge import PivotError, VariantPivoter, merge_samples, _add_mult_alt_flags, create_initial_df, build_pivoter, validate_parameters, rearrange_columns, determine_input_keys, get_headers_and_readers, create_dict, cleanup_df, combine_format_columns, remove_non_jq_tags, add_all_tags, sort_format_tags, determine_merge_execution_context, print_new_execution_context, determine_caller_and_split_mult_alts, validate_samples_for_callers, validate_sample_caller_vcfs, create_new_line, create_merging_dict, remove_old_columns
import jacquard.merge as merge
import jacquard.utils as utils
from jacquard.vcf import FileReader
import jacquard.logger as logger
from jacquard.utils import JQException
import test.mock_module as mock_module
import test_case as test_case

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

mock_log_called = False

def mock_log(msg, *args):
    global mock_log_called
    mock_log_called = True

def dataframe(input_data, sep="\t", index_col=None):
    def tupelizer(thing):
        if isinstance(thing, str) and thing.startswith("(") and thing.endswith(")"):
            return ast.literal_eval(thing)
        return thing

    df = pd.read_csv(StringIO(input_data), sep=sep, header=False, dtype='str',
                     index_col=index_col)
    new_cols = [tupelizer(col) for col in list(df.columns.values)]
    df.columns = pd.core.index.Index(new_cols)

    return df

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
        
class MockReader():
    def __init__(self, content):
        lines = [line + "\n" for line in content.split("\n") if line != ""]
        self._iter = lines.__iter__()
        self.wasClosed = False

    def __iter__(self):
        return self._iter

    def close(self):
        self.wasClosed=True

class MergeTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output
        self.original_info = logger.info
        self.original_error = logger.error
        self.original_warning = logger.warning
        self.original_debug = logger.debug
        self._change_mock_logger()

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr
        self._reset_mock_logger()

    def _change_mock_logger(self):
        global mock_log_called
        mock_log_called = False
        global mock_log
        logger.info = mock_log
        logger.error = mock_log
        logger.warning = mock_log
        logger.debug = mock_log

    def _reset_mock_logger(self):
        logger.info = self.original_info
        logger.error = self.original_error
        logger.warning = self.original_warning
        logger.debug = self.original_debug

    def testExecute_multAltsSplitCorrectly(self):
        vcfRecordFormat = "##jacquard.tag.caller={}\n" + \
            "#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|A|B\n" + \
            "chr1|42|.|{}|{}|.|.|{}|JQ_AF_{}|0.1|0.2\n"
        vcfRecordFormat = vcfRecordFormat.replace("|", "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.mutect.vcf", vcfRecordFormat.format("MuTect","T","A,C","INFO_Mutect","MT"))
            input_dir.write("A.strelka.vcf", vcfRecordFormat.format("Strelka","T","A,C","INFO_Strelka","SK"))
            input_dir.write("A.varscan.vcf", vcfRecordFormat.format("VarScan","T","A,C","INFO_VarScan","VS"))

            args = Namespace(input=input_dir.path, 
                         output=os.path.join(output_dir.path,"tmp.vcf"), 
                         allow_inconsistent_sample_sets=False,
                         keys=None) 
            merge.execute(args, [])

            actual_merged = output_dir.read('tmp.vcf').split("\n")
            
            print actual_merged
            self.assertEquals(14, len(actual_merged))

    def test_addFiles(self):
        rows = ["CHROM", "POS", "REF", "ALT"]

        pivoter = VariantPivoter(rows)
        sample_A_file = \
'''CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSamp1
1\t23\tA\tT\tfoo\tDP:ESAF\t1:0.2
2\t24\tA\tT\tfoo\tDP:ESAF\t12:0.2
3\t25\tA\tT\tfoo\tDP:ESAF\t31:0.2
4\t26\tA\tT\tfoo\tDP:ESAF\t6:0.2'''
        sample_B_file = \
'''CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSamp2
1\t23\tA\tT\tfoo\tDP:ESAF\t5:0.2
2\t24\tA\tT\tfoo\tDP:ESAF\t2:0.2
3\t25\tA\tT\tfoo\tDP:ESAF\t74:0.2
4\t26\tA\tT\tfoo\tDP:ESAF\t25:0.2'''
    
        pivoter.add_file(StringIO(sample_A_file), 0, "MuTect", "file1")
        pivoter.add_file(StringIO(sample_B_file), 0, "MuTect", "file2")
        
        actual_df = pivoter._combined_df
        actual_df.columns.names = [""]

        expected_string = \
'''CHROM\tPOS\tREF\tALT\tINFO\tMuTect|file1|FORMAT\tMuTect|file1|Samp1\tMuTect|file2|FORMAT\tMuTect|file2|Samp2
1\t23\tA\tT\t.\tDP:ESAF\t1:0.2\tDP:ESAF\t5:0.2
2\t24\tA\tT\t.\tDP:ESAF\t12:0.2\tDP:ESAF\t2:0.2
3\t25\tA\tT\t.\tDP:ESAF\t31:0.2\tDP:ESAF\t74:0.2
4\t26\tA\tT\t.\tDP:ESAF\t6:0.2\tDP:ESAF\t25:0.2'''
        expected_df = dataframe(expected_string)
        expected_df.columns.names = [""]

        tm.assert_frame_equal(expected_df, actual_df)
         
    def test_validateSampleCallerVcfs(self):
        dataString1 = \
'''COORDINATE\tVarScan|foo|FORMAT\tVarScan|sample_A\tVarScan|sample_B\tMuTect|foo|FORMAT\tMuTect|foo|sample_A\tMuTect|foo|sample_A
1\tGT:ESAF\t10:0.2\t100:0.2
2\tGT:ESAF\t20:0.2\t200:0.2
3\tGT:ESAF\t30:0.2\t300:0.2
4\tGT:ESAF\t40:0.2\t400:0.2'''
        df = pd.read_csv(StringIO(dataString1), sep="\t", header=False, dtype='str', mangle_dupe_cols=False)

        self.assertRaisesRegexp(JQException,
                                    "Some samples have calls for the same caller in more than one file. Adjust or move problem input files and try again.",
                                    merge.validate_sample_caller_vcfs,
                                    df)
#         self.assertTrue("Sample [foo|sample_A] appears to be called by [MuTect] in multiple files." in self.output.getvalue())
        global mock_log_called
        self.assertTrue(mock_log_called)
        
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
        
    def Ytest_validateSampleData_nonUniqueRows(self):
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
        
    #TODO: move out of reference files
    def test_determineInputKeysVcf(self):
        input_dir = TEST_DIRECTORY + "/functional_tests/test_input/test_input_keys_vcf"
        actual_lst = determine_input_keys(input_dir)
        
        expected_lst = ["CHROM", "POS", "REF", "ALT"]
        
        self.assertEquals(expected_lst, actual_lst)
        
    def test_determineInputKeysInvalid(self):
        input_dir = TEST_DIRECTORY + "/functional_tests/test_input/test_input_keys_invalid"
        
        self.assertRaises(PivotError, determine_input_keys, input_dir)
    
    ##get headers, readers
    def test_getHeadersAndReaders(self):
        input_dir = TEST_DIRECTORY + "/functional_tests/test_input/test_input_valid"
        in_files = sorted(glob.glob(os.path.join(input_dir,"*")))
        sample_file_readers, headers, header_names, first_line, meta_headers = get_headers_and_readers(in_files)
        
        self.assertEquals([os.path.join(input_dir, "foo1.txt")], sample_file_readers)
        self.assertEquals([3], headers)
        self.assertEquals("CHROM\tPOS\tREF\tALT\tGENE_SYMBOL\tFORMAT\tSample_2384\tSample_2385", header_names.rstrip())
        self.assertEquals("1\t2342\tA\tT\tEGFR\tGT:DP\t1/1:241\t0/1:70", first_line[0].rstrip())
        self.assertEquals(['##FORMAT=<ID=JQ_FOO'], meta_headers)
        
    def test_getHeadersAndReaders_invalid(self):
        input_dir = TEST_DIRECTORY + "/functional_tests/test_input/test_input_keys_txt"
        in_files = sorted(glob.glob(os.path.join(input_dir,"*")))
        
        self.assertRaisesRegexp(JQException,
                                r"VCF file\(s\) .* have no Jacquard tags. Run \[jacquard tag\] on these files and try again.",
                                get_headers_and_readers,
                                in_files)
                
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
'''#CHROM\tPOS\tREF\tINFO
1\t42\tA\t.
2\t43\tC\t.''');

        expected_df = pd.read_csv(expected_data, sep="\t", header=False, dtype='str')

        tm.assert_frame_equal(expected_df, actual_df)

    def test_mergeSamples_emptyCombinedDf(self): 
        dataString = \
'''CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tsample_A\tsample_B
1\t23\tA\tG\tfoo\tGT:ESAF\t10:0.2\t100:0.2
2\t24\tA\tG\tfoo\tGT:ESAF\t20:0.2\t200:0.2
3\t25\tA\tG\tfoo\tGT:ESAF\t30:0.2\t300:0.2
4\t26\tA\tG\tfoo\tGT:ESAF\t40:0.2\t400:0.2'''
        df = pd.read_csv(StringIO(dataString), sep="\t", header=False, dtype='str')

        combined_df = pd.DataFrame()
        actual_df = merge_samples(df, combined_df, ["CHROM", "POS", "REF", "ALT"])

        expected_data = StringIO(
'''CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tsample_A\tsample_B
1\t23\tA\tG\tfoo\tGT:ESAF\t10:0.2\t100:0.2
2\t24\tA\tG\tfoo\tGT:ESAF\t20:0.2\t200:0.2
3\t25\tA\tG\tfoo\tGT:ESAF\t30:0.2\t300:0.2
4\t26\tA\tG\tfoo\tGT:ESAF\t40:0.2\t400:0.2''')
        expected_df = pd.read_csv(expected_data, sep="\t", header=False, dtype='str')
        
        tm.assert_frame_equal(expected_df, actual_df)
    
    def test_addMultAltFlags(self):
        dataString = \
'''CHROM\tPOS\tREF\tALT\tINFO\tfile1|FORMAT\tfile1|sample_A
1\t12\tA\tG\t.\tGT:ESAF\t10:0.2
1\t12\tA\tC\t.\tGT:ESAF\t10:0.2'''
        df = pd.read_csv(StringIO(dataString), sep="\t", header=False, dtype='str')
        
        actual_df = _add_mult_alt_flags(df)
        
        expectedString =  \
'''CHROM\tPOS\tREF\tALT\tINFO\tfile1|FORMAT\tfile1|sample_A
1\t12\tA\tG\tMult_Alt\tGT:ESAF\t10:0.2
1\t12\tA\tC\tMult_Alt\tGT:ESAF\t10:0.2'''
        expected_df = pd.read_csv(StringIO(expectedString), sep="\t", header=False, dtype='str')
        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_mergeSamples_multAlts(self):
        dataString1 = \
'''CHROM\tPOS\tREF\tALT\tINFO\tfile1|FORMAT\tfile1|sample_A\tfile1|sample_B
1\t12\tA\tG\t.\tGT:ESAF\t10:0.2\t100:0.2'''
        df1 = pd.read_csv(StringIO(dataString1), sep="\t", header=False, dtype='str')
        
        dataString2 = \
'''CHROM\tPOS\tREF\tALT\tINFO\tfile2|FORMAT\tfile2|sample_C\tfile2|sample_D
1\t12\tA\tC\t.\tGT:ESAF\t10:0.2\t100:0.2'''
        df2 = pd.read_csv(StringIO(dataString2), sep="\t", header=False, dtype='str')

        combined_df = df2
        actual_df = merge_samples(combined_df, df1, ["CHROM", "POS", "REF", "ALT"])

        expectedString = \
'''CHROM\tPOS\tREF\tALT\tINFO\tfile1|FORMAT\tfile1|sample_A\tfile1|sample_B\tfile2|FORMAT\tfile2|sample_C\tfile2|sample_D
1\t12\tA\tG\t.\tGT:ESAF\t10:0.2\t100:0.2\tnan\tnan\tnan
1\t12\tA\tC\t.\tnan\tnan\tnan\tGT:ESAF\t10:0.2\t100:0.2'''
        expected_df = pd.read_csv(StringIO(expectedString), sep="\t", header=False, dtype='str')
        
        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_mergeSamples_populatedCombinedDf(self): 
        dataString1 = \
'''CHROM\tPOS\tREF\tALT\tINFO\tfile1|FORMAT\tfile1|sample_A\tfile1|sample_B
1\t12\tA\tG\tfoo\tGT:ESAF\t10:0.2\t100:0.2
2\t13\tA\tG\tfoo\tGT:ESAF\t20:0.2\t200:0.2
3\t14\tA\tG\tfoo\tGT:ESAF\t30:0.2\t300:0.2
4\t15\tA\tG\tfoo\tGT:ESAF\t40:0.2\t400:0.2'''
        df1 = pd.read_csv(StringIO(dataString1), sep="\t", header=False, dtype='str')
        
        dataString2 = \
'''CHROM\tPOS\tREF\tALT\tINFO\tfile2|FORMAT\tfile2|sample_C\tfile2|sample_D
1\t12\tA\tG\tfoo\tGT:ESAF\t10:0.2\t100:0.2
2\t13\tA\tG\tfoo\tGT:ESAF\t20:0.2\t200:0.2
3\t14\tA\tG\tfoo\tGT:ESAF\t30:0.2\t300:0.2
4\t15\tA\tG\tfoo\tGT:ESAF\t40:0.2\t400:0.2'''
        df2 = pd.read_csv(StringIO(dataString2), sep="\t", header=False, dtype='str')

        combined_df = df2
        actual_df = merge_samples(combined_df, df1, ["CHROM", "POS", "REF", "ALT"])

        expected_data = StringIO(
'''CHROM\tPOS\tREF\tALT\tINFO\tfile1|FORMAT\tfile1|sample_A\tfile1|sample_B\tfile2|FORMAT\tfile2|sample_C\tfile2|sample_D
1\t12\tA\tG\tfoo\tGT:ESAF\t10:0.2\t100:0.2\tGT:ESAF\t10:0.2\t100:0.2
2\t13\tA\tG\tfoo\tGT:ESAF\t20:0.2\t200:0.2\tGT:ESAF\t20:0.2\t200:0.2
3\t14\tA\tG\tfoo\tGT:ESAF\t30:0.2\t300:0.2\tGT:ESAF\t30:0.2\t300:0.2
4\t15\tA\tG\tfoo\tGT:ESAF\t40:0.2\t400:0.2\tGT:ESAF\t40:0.2\t400:0.2''')
        expected_df = pd.read_csv(expected_data, sep="\t", header=False, dtype='str')
        
        print actual_df
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
    
    def test_determineCaller_valid(self):
        reader = MockReader("##foo\n##jacquard.tag.caller=MuTect\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n1\t2324\t.\tA\tG,T\t.\t.\t.\tJQ_AF_VS:DP:FOO\t0.234,0.124:78:25,312")
        writer = MockWriter()
        unknown_callers = 0
        caller, unknown_callers = determine_caller_and_split_mult_alts(reader, writer, unknown_callers)
        self.assertEquals("MuTect", caller)
        self.assertEquals(0, unknown_callers)

        self.assertEquals(["##foo", "##jacquard.tag.caller=MuTect", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2","1\t2324\t.\tA\tG\t.\t.\t.\tJQ_AF_VS:DP:FOO\t0.234:78:25,312","1\t2324\t.\tA\tT\t.\t.\t.\tJQ_AF_VS:DP:FOO\t0.124:78:25,312"], writer.lines())
     
    def test_determineCaller_invalid(self):
        reader = MockReader("##foo\n##jacquard.tag.foo\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n1\t2324\t.\tA\tG,T\t.\t.\t.\tJQ_AF_VS:DP:FOO\t0.234,0.124:78:25,312")
        writer = MockWriter()
        unknown_callers = 2
        caller, unknown_callers = determine_caller_and_split_mult_alts(reader, writer, unknown_callers)
        self.assertEquals("unknown", caller)
        self.assertEquals(3, unknown_callers)
        self.assertEquals(["##foo", "##jacquard.tag.foo", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2","1\t2324\t.\tA\tG\t.\t.\t.\tJQ_AF_VS:DP:FOO\t0.234:78:25,312","1\t2324\t.\tA\tT\t.\t.\t.\tJQ_AF_VS:DP:FOO\t0.124:78:25,312"], writer.lines())
        
    def test_createNewLine(self):
        fields = ["1", "42", ".", "A", "G,CT", ".", ".", ".", "DP:JQ_VS_AF:AF", "23:0.24,0.32:0.2354,0.324", "23:0.25,0.36:0.254,0.3456"]

        alt_allele_number = 0
        actual_line = create_new_line(alt_allele_number, fields)
        expected_line = "\t".join(["1", "42", ".", "A", "G", ".", ".", ".", "DP:JQ_VS_AF:AF", "23:0.24:0.2354,0.324", "23:0.25:0.254,0.3456\n"])
        self.assertEquals(expected_line, actual_line)

        alt_allele_number = 1
        actual_line = create_new_line(alt_allele_number, fields)
        expected_line = "\t".join(["1", "42", ".", "A", "CT", ".", ".", ".", "DP:JQ_VS_AF:AF", "23:0.32:0.2354,0.324", "23:0.36:0.254,0.3456\n"])
        self.assertEquals(expected_line, actual_line)

    def test_validateSamplesForCallers_valid(self):
        context = ["jacquard.foo=file1|13523|tumor(file1.vcf)", "jacquard.foo=file2|13523|tumor(file2.vcf)", "jacquard.foo=file3|13523|tumor(file3.vcf)"]
        value = validate_samples_for_callers(context, False)
        self.assertEqual(1, value)
        
    def test_validateSamplesForCallers_validAllow(self):
        context = ["jacquard.foo=file1|13523|tumor(file1.vcf)", "jacquard.foo=file2|13523|tumor(file2.vcf)", "jacquard.foo=file3|13523|tumor(file3.vcf)"]
        value = validate_samples_for_callers(context, True)
        self.assertEqual(1, value)
    
    def test_validateSamplesForCallers_invalid(self):
        context = ["jacquard.foo=file1|13523|tumor(file1.vcf)", "jacquard.foo=file2|23|tumor(file2.vcf)", "jacquard.foo=file3|768|tumor(file3.vcf)"]
    
        self.assertRaisesRegexp(JQException,
                                "Some samples were not present for all callers. Review log warnings and move/adjust input files as appropriate.",
                                validate_samples_for_callers,
                                context, False)
        
    def test_validateSamplesForCallers_invalidAllow(self):
        context = ["jacquard.foo=file1|13523|tumor(file1.vcf)", "jacquard.foo=file2|23|tumor(file2.vcf)", "jacquard.foo=file3|768|tumor(file3.vcf)"]
        value = validate_samples_for_callers(context, True)
        self.assertEqual(1, value)
        
    def test_createDict(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tMuTect|file1|FORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tDP\t57\t57
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tDP\t58\t57
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tDP\t59\t57
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tDP\t60\t57
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tDP\t60\t57
'''
        df = dataframe(input_string)
        row = 0
        col = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "MuTect|file1|FORMAT", "MuTect|file1|sample_A", "MuTect|file1|sample_B"]
        expected_file_dict = {'MuTect': [OrderedDict([('DP', '57'), ('sample_name', 'MuTect|file1|sample_A')]), OrderedDict([('DP', '57'), ('sample_name', 'MuTect|file1|sample_B')])]}
        file_dict, all_tags = create_dict(df, row, col)
        self.assertEqual(expected_file_dict, file_dict)
    
    def test_createDict_multFiles(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tMuTect|file1|FORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B\tMuTect|file2|FORMAT\tMuTect|file2|sample_A\tMuTect|file2|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tDP\t57\t57\tDP\t57\t57
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tDP\t58\t57\tDP\t57\t57
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tDP\t59\t57\tDP\t57\t57
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tDP\t60\t57\tDP\t57\t57
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tDP\t60\t57\tDP\t57\t57
'''
        df = dataframe(input_string)
        row = 0
        col = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "MuTect|file1|FORMAT", "MuTect|file1|sample_A", "MuTect|file1|sample_B", "MuTect|file2|FORMAT", "MuTect|file2|sample_A", "MuTect|file2|sample_B"]
        expected_file_dict = {'MuTect': [OrderedDict([('DP', '57'), ('sample_name', 'MuTect|file1|sample_A')]), OrderedDict([('DP', '57'), ('sample_name', 'MuTect|file1|sample_B')]), OrderedDict([('DP', '57'), ('sample_name', 'MuTect|file2|sample_A')]), OrderedDict([('DP', '57'), ('sample_name', 'MuTect|file2|sample_B')])]}
        file_dict, all_tags = create_dict(df, row, col)
        self.assertEqual(expected_file_dict, file_dict)
        
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
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B\tMuTect|file2|FORMAT\tMuTect|file2|sample_A\tMuTect|file2|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tsample_name:DP\t57:MuTect|file1|sample_A\t57:MuTect|file1|sample_B\tDP:sample_name\t57:MuTect|file2|sample_A\t57:MuTect|file2|sample_B
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tsample_name:DP\t58:MuTect|file1|sample_A\t57:MuTect|file1|sample_B\tDP:sample_name\t57:MuTect|file2|sample_A\t57:MuTect|file2|sample_B
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tsample_name:DP\t59:MuTect|file1|sample_A\t57:MuTect|file1|sample_B\tDP:sample_name\t57:MuTect|file2|sample_A\t57:MuTect|file2|sample_B
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tsample_name:DP\t60:MuTect|file1|sample_A\t57:MuTect|file1|sample_B\tDP:sample_name\t57:MuTect|file2|sample_A\t57:MuTect|file2|sample_B
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tsample_name:nan\tnan:MuTect|file1|sample_A\tnan:MuTect|file1|sample_B\tDP:sample_name\t57:MuTect|file2|sample_A\t57:MuTect|file2|sample_B
'''
        df = dataframe(input_string)
        file_dict = {"MuTect|file1":"foo", "MuTect|file2": "foo"}
        actual_df = cleanup_df(df, file_dict)
        
        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B\tMuTect|file2|sample_A\tMuTect|file2|sample_B
1\t2\t.\tA\tG\t.\t.\tfoo\tDP\t57\t57\t57\t57
1\t3\t.\tC\tG\t.\t.\tfoo\tDP\t58\t57\t57\t57
8\t4\t.\tA\tT\t.\t.\tfoo\tDP\t59\t57\t57\t57
13\t5\t.\tT\tAA\t.\t.\tfoo\tDP\t60\t57\t57\t57
13\t5\t.\tT\tAAAA\t.\t.\tfoo\t.\t.\t.\t57\t57
'''
        expected_df = dataframe(expected_string)
        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_createMergingDict(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tJQ_DP\t57\t57
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tJQ_DP\t58\t57
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tJQ_DP\t59\t57
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tJQ_DP\t60\t57
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tnan\tnan\tnan
'''
        df = dataframe(input_string)
        row = 0
        col = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "MuTect|file1|sample_A", "MuTect|file1|sample_B"]
        file_dict = create_merging_dict(df, row, col)
        expected_file_dict = {'file1|sample_A': [OrderedDict([('JQ_DP', '57')])], 'file1|sample_B': [OrderedDict([('JQ_DP', '57')])]}
        self.assertEquals(expected_file_dict, file_dict)
    
    def test_createMergingDict_multipleCallers(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B\tVarScan|file1|sample_A\tVarScan|file1|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tJQ_DP_VS:JQ_DP_MT\t.:57\t.:57\t57:.\t57:.
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\JQ_DP_VS:JQ_DP_MT\t.:58\t.:57\t57:.\t57:.
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\JQ_DP_VS:JQ_DP_MT\t.:59\t.:57\t57:.\t57:.
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\JQ_DP_VS:JQ_DP_MT\t.:60\t.:57\t57:.\t57:.
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tnan:JQ_DP_MT\t.:nan\t.:nan\t57:.\t57:.
'''
        df = dataframe(input_string)
        row = 0
        col = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "MuTect|file1|sample_A", "MuTect|file1|sample_B", "VarScan|file1|sample_A", "VarScan|file1|sample_B"]
        file_dict = create_merging_dict(df, row, col)
        expected_file_dict = {'file1|sample_A': [OrderedDict([('JQ_DP_VS', '.'), ('JQ_DP_MT', '57')]), OrderedDict([('JQ_DP_VS', '57'), ('JQ_DP_MT', '.')])], 'file1|sample_B': [OrderedDict([('JQ_DP_VS', '.'), ('JQ_DP_MT', '57')]), OrderedDict([('JQ_DP_VS', '57'), ('JQ_DP_MT', '.')])]}
        self.assertEquals(expected_file_dict, file_dict)
    
    def test_removeOldColumns(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t57\t57\t57\t57
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t58\t57\t57\t57
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t59\t57\t57\t57
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t60\t57\t57\t57
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tJQ_DP_MT\tnan\tnan\t57\t57
'''
        df = dataframe(input_string)
        actual_df = remove_old_columns(df)
        
        expected_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t57\t57
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t57\t57
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t57\t57
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t57\t57
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t57\t57
'''
        expected_df = dataframe(expected_string)
        
        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_combineFormatColumns(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tMuTect|file1|FORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tJQ_DP:JQ_AF\t57:0.23\t57:0.4
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tJQ_DP\t58\t57
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tJQ_DP\t59\t57
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tJQ_DP\t60\t57
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tnan\tnan\tnan
'''
        df = dataframe(input_string)

        actual_df = merge.combine_format_columns(df, 0)
        
        global mock_log_called
        self.assertTrue(mock_log_called)

        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\t.\t.\tfoo\tJQ_AF:JQ_DP\t0.23:57\t0.4:57
1\t3\t.\tC\tG\t.\t.\tfoo\tJQ_DP\t58\t57
8\t4\t.\tA\tT\t.\t.\tfoo\tJQ_DP\t59\t57
13\t5\t.\tT\tAA\t.\t.\tfoo\tJQ_DP\t60\t57
13\t5\t.\tT\tAAAA\t.\t.\tfoo\t.\t.\t.
'''
        expected_df = dataframe(expected_string)

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_combineFormatColumns_differentSampleNames(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tMuTect|file1|FORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B\tVarScan|file1|FORMAT\tVarScan|file1|sample_a\tVarScan|file1|sample_b
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t1\t20\tJQ_DP_VS\t6\t25
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t2\t21\tJQ_DP_VS\t7\t26
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t3\t22\tJQ_DP_VS\t8\t27
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t4\t23\tJQ_DP_VS\t9\t28
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t5\t24\tJQ_DP_VS\t10\t29
'''
        df = dataframe(input_string)


        actual_df = merge.combine_format_columns(df, 0)
        
        global mock_log_called
        self.assertTrue(mock_log_called)
        
        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile1|sample_b\tfile1|sample_A\tfile1|sample_B\tfile1|sample_a
1\t2\t.\tA\tG\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t.:25\t1:.\t20:.\t.:6
1\t3\t.\tC\tG\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t.:26\t2:.\t21:.\t.:7
8\t4\t.\tA\tT\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t.:27\t3:.\t22:.\t.:8
13\t5\t.\tT\tAA\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t.:28\t4:.\t23:.\t.:9
13\t5\t.\tT\tAAAA\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t.:29\t5:.\t24:.\t.:10
'''
        expected_df = dataframe(expected_string, "\t", False)

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_combineFormatColumns_missingTags(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tMuTect|file1|FORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B\tVarScan|file1|FORMAT\tVarScan|file1|sample_A\tVarScan|file1|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t1\t20\t.\t.\t.
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t2\t21\tJQ_DP_VS\t7\t26
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\t.\t.\t.\tJQ_DP_VS\t8\t27
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t4\t23\tJQ_DP_VS\t9\t28
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t5\t24\tJQ_DP_VS\t10\t29
'''
        
        df = dataframe(input_string)

        actual_df = merge.combine_format_columns(df, 0)
        
        global mock_log_called
        self.assertTrue(mock_log_called)
        
        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\t.\t.\tfoo\tJQ_DP_MT\t1\t20
1\t3\t.\tC\tG\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t2:7\t21:26
8\t4\t.\tA\tT\t.\t.\tfoo\tJQ_DP_VS\t8\t27
13\t5\t.\tT\tAA\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t4:9\t23:28
13\t5\t.\tT\tAAAA\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t5:10\t24:29
'''

        expected_df = dataframe(expected_string)

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_combineFormatColumns_validateOrder(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tMuTect|file1|FORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B\tVarScan|file1|FORMAT\tVarScan|file1|sample_A\tVarScan|file1|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t1\t20\tJQ_DP_VS\t6\t25
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t2\t21\tJQ_DP_VS\t7\t26
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t3\t22\tJQ_DP_VS\t8\t27
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t4\t23\tJQ_DP_VS\t9\t28
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tJQ_DP_MT\t5\t24\tJQ_DP_VS\t10\t29
'''
        df = dataframe(input_string)

        actual_df = merge.combine_format_columns(df, 0)
        
        global mock_log_called
        self.assertTrue(mock_log_called)
        
        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t1:6\t20:25
1\t3\t.\tC\tG\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t2:7\t21:26
8\t4\t.\tA\tT\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t3:8\t22:27
13\t5\t.\tT\tAA\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t4:9\t23:28
13\t5\t.\tT\tAAAA\t.\t.\tfoo\tJQ_DP_MT:JQ_DP_VS\t5:10\t24:29
'''
        expected_df = dataframe(expected_string)

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_combineFormatColumns_inconsistentSampleSets(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tMuTect|file1|FORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B\tMuTect|file2|FORMAT\tMuTect|file2|sample_A\tMuTect|file2|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t57:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t58:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t59:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t60:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tnan\tnan\tnan\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
'''
        df = dataframe(input_string)

        actual_df = merge.combine_format_columns(df, 1)
        
        global mock_log_called
        self.assertTrue(mock_log_called)

        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile2|sample_A\tfile2|sample_B\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1
1\t3\t.\tC\tG\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1
8\t4\t.\tA\tT\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1
13\t5\t.\tT\tAA\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1
13\t5\t.\tT\tAAAA\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t.\t.
'''
        
        expected_df = dataframe(expected_string)

        tm.assert_frame_equal(expected_df, actual_df)

    def test_combineFormatColumns_consistentSampleSets(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tMuTect|file1|FORMAT\tMuTect|file1|sample_A\tMuTect|file1|sample_B\tMuTect|file2|FORMAT\tMuTect|file2|sample_A\tMuTect|file2|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t57:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
1\t3\t.\tC\tG\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t58:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
8\t4\t.\tA\tT\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t59:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
13\t5\t.\tT\tAA\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t60:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
13\t5\t.\tT\tAAAA\tQUAL\tFILTER\tfoo\tDP:JQ_GT\t61:0/1\t57:0/1\tDP:AF:JQ_GT\t57:0.2:0/1\t57:0.2:0/1
'''
        df = dataframe(input_string)


        actual_df = merge.combine_format_columns(df, 0)
        
        global mock_log_called
        self.assertTrue(mock_log_called)

        expected_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile2|sample_A\tfile2|sample_B\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1
1\t3\t.\tC\tG\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1
8\t4\t.\tA\tT\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1
13\t5\t.\tT\tAA\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1
13\t5\t.\tT\tAAAA\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1
'''
        
        expected_df = dataframe(expected_string)

        tm.assert_frame_equal(expected_df, actual_df)

    def test_combineFormatColumns_differingVariants(self):
        input_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tMuTect|file1|FORMAT\tMuTect|file1|sample_A\tMuTect|file2|FORMAT\tMuTect|file2|sample_B
1\t2\t.\tA\tG\tQUAL\tFILTER\tfoo\tnan\tnan\tJQ_DP:JQ_GT\t57:0/1
1\t3\t.\tA\tG\tQUAL\tFILTER\tfoo\tJQ_DP:JQ_GT\t58:0/1\tnan\tnan
'''
        df = dataframe(input_string)

        actual_df = merge.combine_format_columns(df, 0)
        
        global mock_log_called
        self.assertTrue(mock_log_called)

        expected_string = \
'''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile2|sample_B\tfile1|sample_A
1\t2\t.\tA\tG\t.\t.\tfoo\tJQ_DP:JQ_GT\t57:0/1\t.:.
1\t3\t.\tA\tG\t.\t.\tfoo\tJQ_DP:JQ_GT\t.:.\t58:0/1
'''
        
        expected_df = dataframe(expected_string)

        tm.assert_frame_equal(expected_df, actual_df)
        
    def test_removeNonJQTags(self):
        fake_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile2|sample_A\tfile2|sample_B\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1'''
        df = dataframe(fake_string)
        file_dict = {'file2': [OrderedDict([('DP', '57'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('DP', '57'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_B')])]}
        actual_file_dict = remove_non_jq_tags(df, file_dict)
        expected_file_dict = {'file2': [OrderedDict([('JQ_foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('JQ_foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('JQ_foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('JQ_foo', '.'), ('sample_name', 'file1|sample_B')])]}
       
        self.assertEquals(expected_file_dict, actual_file_dict)
    
    def test_removeNonJQTags_noTags(self):
        fake_string = \
        '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfile2|sample_A\tfile2|sample_B\tfile1|sample_A\tfile1|sample_B
1\t2\t.\tA\tG\t.\t.\tfoo\tJQ_GT\t0/1\t0/1\t0/1\t0/1'''
        df = dataframe(fake_string)
        file_dict = {'file2': [OrderedDict([('DP', '57'), ('foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('DP', '57'), ('foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('DP', '57'), ('foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('DP', '57'), ('foo', '.'), ('sample_name', 'file1|sample_B')])]}
        actual_file_dict = remove_non_jq_tags(df, file_dict)
        expected_file_dict = {'file2': [OrderedDict([('sample_name', 'file2|sample_A')]), OrderedDict([('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('sample_name', 'file1|sample_A')]), OrderedDict([('sample_name', 'file1|sample_B')])]}
       
        self.assertEquals(expected_file_dict, actual_file_dict)
    
    def test_addAllTags(self):
        file_dict = {'file2': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('JQ_foo', '1'), ('sample_name', 'file2|sample_B')])], 'file1': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_B')])]}
        sample_keys = ["JQ_DP", "JQ_foo"]
        
        actual_file_dict = add_all_tags(file_dict, sample_keys)
        expected_file_dict = {'file2': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '1'), ('sample_name', 'file2|sample_A')]), OrderedDict([('JQ_foo', '1'), ('sample_name', 'file2|sample_B'), ('JQ_DP', '^')])], 'file1': [OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_A')]), OrderedDict([('JQ_DP', '57'), ('JQ_foo', '.'), ('sample_name', 'file1|sample_B')])]}
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
        sample_file = "file1.vcf"
        count = 1
        all_merge_context, all_merge_column_context = determine_merge_execution_context(all_merge_context, all_merge_column_context, sample_columns, sample_file, count)

        self.assertEqual(["##jacquard.merge.file1=file1.vcf(['samp1', 'samp2'])"], all_merge_context)
        self.assertEqual(['##jacquard.merge.sample_column1=file1|samp1(file1.vcf)', '##jacquard.merge.sample_column2=file1|samp2(file1.vcf)'], all_merge_column_context)

    def test_printNewExecutionContext(self):
        out_file = MockWriter()
        execution_context = ["##jacquard", "##foo_bar", "##baz"]
        print_new_execution_context(execution_context, out_file)
        lines = out_file.lines()
        self.assertEqual(['##jacquard', '##foo_bar', '##baz'], lines)
        self.assertEqual(True, out_file.wasClosed)

class MergeFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_merge(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "04_merge")
            input_dir = os.path.join(module_testdir, "input")
            output_file = os.path.join(output_dir.path, "merged.vcf")

            command = ["merge", input_dir, output_file, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

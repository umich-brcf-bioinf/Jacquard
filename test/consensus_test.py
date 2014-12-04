# pylint: disable=C0111
from collections import OrderedDict
import numpy
from StringIO import StringIO
import sys
import unittest
from argparse import Namespace
import jacquard.utils as utils
import glob
from testfixtures import TempDirectory
import os
from jacquard.consensus import iterate_file, add_consensus, process_line, calculate_consensus, create_consensus_dict, get_consensus_som, get_consensus, add_zscore, calculate_zscore
import jacquard.consensus as consensus
import jacquard.logger as logger
from jacquard.vcf import FileReader

mock_log_called = False

def mock_log(msg, *args):
    global mock_log_called
    mock_log_called = True

class ConsensusTestCase(unittest.TestCase):
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

    def test_iterateFile_consensus(self):
        reader = MockReader("##FOO\n#CHROM\n1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_VS_HC_SOM:DP:JQ_VS_AF:{0}AF:{0}DP\t1:0:0.0:0.0:0\t0:1:1.2:1.2:0\t1:3:1.0:1.0:0\n".format(consensus.JQ_CONSENSUS_TAG))
        af_range = []
        dp_range = []
        function = "consensus"
        meta_headers, header, lines = iterate_file(reader, af_range, dp_range, function)

        self.assertEquals(["##FOO\n"], meta_headers)
        self.assertEquals("#CHROM\n", header)
        self.assertEquals(["1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_VS_HC_SOM:DP:JQ_VS_AF:{0}AF:{0}DP:{0}SOM_SUM:{0}AF_RANGE:{0}DP_RANGE\t1:0:0.0:0.0:0.0:1:0.0:0.0\t0:1:1.2:1.2:0.0:0:0.0:0.0\t1:3:1.0:1.0:0.0:1:0.0:0.0\n".format(consensus.JQ_CONSENSUS_TAG)], lines)

    def test_iterateFile_zscore(self):
        reader = MockReader("##FOO\n#CHROM\n1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_HC_SOM_VS:DP:JQ_AF_VS:{0}AF:{0}DP:{0}SOM_SUM:{0}AF_RANGE:{0}DP_RANGE\t1:0:0.0:0.0:0.0:1:0.0:0\t0:1:1.2:1.2:0.0:0:0.0:0\t1:3:1.0:1.0:0.0:1:0.0:0\n".format(consensus.JQ_CONSENSUS_TAG))
        af_range = [0.2, 0.45]
        dp_range = [45.1, 86.0]
        function = "zscore"
        meta_headers, header, lines = iterate_file(reader, af_range, dp_range, function)

        self.assertEquals(["##FOO\n"], meta_headers)
        self.assertEquals("#CHROM\n", header)
        self.assertEquals(["1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_HC_SOM_VS:DP:JQ_AF_VS:{0}AF:{0}DP:{0}SOM_SUM:{0}AF_RANGE:{0}DP_RANGE:{0}AF_RANGE_ZSCORE:{0}DP_RANGE_ZSCORE\t1:0:0.0:0.0:0.0:1:0.0:0:-2.6:-3.21\t0:1:1.2:1.2:0.0:0:0.0:0:-2.6:-3.21\t1:3:1.0:1.0:0.0:1:0.0:0:-2.6:-3.21\n".format(consensus.JQ_CONSENSUS_TAG)], lines)

    def test_addConsensusSomatic(self):
        meta_headers = []
        header = "#CHROM"
        lines = ["1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_AF_VS:{0}AF:{0}SOM_SUM:{0}DP:{0}AF_RANGE:{0}DP_RANGE\t1:0:0.0:0.0:1:0:0:0\t0:1:1.2:1.2:0:0:0:0\t1:3:1.0:1.0:1:0:0:0".format(consensus.JQ_CONSENSUS_TAG)]
        writer = MockWriter()
        output_file = "output"
        add_consensus(meta_headers, header, lines, writer, output_file)
        self.assertEquals(['##FORMAT=<ID={0}SOM_SUM,Number=1,Type=Integer,Description="Jacquard consensus somatic call = sum(*HC_SOM*)">'.format(consensus.JQ_CONSENSUS_TAG),
                            '##FORMAT=<ID={0}AF,Number=A,Type=Integer,Description="Jacquard consensus somatic call = average(*AF*)">'.format(consensus.JQ_CONSENSUS_TAG),
                            '##FORMAT=<ID={0}DP,Number=1,Type=Integer,Description="Jacquard consensus depth = average(*DP*)">'.format(consensus.JQ_CONSENSUS_TAG),
                            '#CHROM',
                            '1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_AF_VS:{0}AF:{0}SOM_SUM:{0}DP:{0}AF_RANGE:{0}DP_RANGE\t1:0:0.0:0.0:1:0:0:0\t0:1:1.2:1.2:0:0:0:0\t1:3:1.0:1.0:1:0:0:0'.format(consensus.JQ_CONSENSUS_TAG)], writer.lines())

    def test_addConsensusSomatic_multipleCallers(self):
        meta_headers = []
        header = "#CHROM"
        lines = ["1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_SOM_MT:{0}AF:{0}SOM_SUM:{0}DP:{0}AF_RANGE:{0}DP_RANGE\t1:0:1:0:2:0:0:0\t0:1:0:0:0:0:0:0\t1:3:0:0:1:0:0:0\n".format(consensus.JQ_CONSENSUS_TAG)]
        writer = MockWriter()
        output_file = "output"
        add_consensus(meta_headers, header, lines, writer, output_file)
        self.assertEquals(['##FORMAT=<ID={0}SOM_SUM,Number=1,Type=Integer,Description="Jacquard consensus somatic call = sum(*HC_SOM*)">'.format(consensus.JQ_CONSENSUS_TAG),
                           '##FORMAT=<ID={0}AF,Number=A,Type=Integer,Description="Jacquard consensus somatic call = average(*AF*)">'.format(consensus.JQ_CONSENSUS_TAG),
                           '##FORMAT=<ID={0}DP,Number=1,Type=Integer,Description="Jacquard consensus depth = average(*DP*)">'.format(consensus.JQ_CONSENSUS_TAG),
                           '#CHROM', 
                           '1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_SOM_MT:{0}AF:{0}SOM_SUM:{0}DP:{0}AF_RANGE:{0}DP_RANGE\t1:0:1:0:2:0:0:0\t0:1:0:0:0:0:0:0\t1:3:0:0:1:0:0:0'.format(consensus.JQ_CONSENSUS_TAG)], writer.lines())

    def test_calculateConsensus(self):
        combined_dict = OrderedDict([("DP", "34"), ("JQ_VS_HC_SOM", "1"), ("JQ_MT_HC_SOM", "0"), ("JQ_VS_AF", "0.5"), ("JQ_VS_DP", "53")])
        actual_dict, af, dp = calculate_consensus(combined_dict, [], [])
        expected_dict = OrderedDict([("DP", "34"), ("JQ_VS_HC_SOM", "1"), ("JQ_MT_HC_SOM", "0"), ("JQ_VS_AF", "0.5"), ("JQ_VS_DP", "53"), ("{0}AF".format(consensus.JQ_CONSENSUS_TAG), "0.5"), ("{0}SOM_SUM".format(consensus.JQ_CONSENSUS_TAG), "1"), ("{0}DP".format(consensus.JQ_CONSENSUS_TAG), "53.0"), ('{0}AF_RANGE'.format(consensus.JQ_CONSENSUS_TAG), '0.0'), ('{0}DP_RANGE'.format(consensus.JQ_CONSENSUS_TAG), '0.0')])
        self.assertEquals(expected_dict, actual_dict)
        self.assertEquals([0], af)
        self.assertEquals([0], dp)

    def test_calculateConsensus_multipleCallers(self):
        combined_dict = OrderedDict([("DP", "34"), ("JQ_VS_HC_SOM", "1"), ("JQ_MT_HC_SOM", "1"), ("JQ_FOO_HC_SOM", "1"), ("JQ_VS_AF", "0.5"), ("JQ_MT_AF", "0.43"), ("JQ_MT_DP", "32"), ("JQ_VS_DP", "23")])
        actual_dict, af, dp = calculate_consensus(combined_dict, [], [])
        expected_dict = OrderedDict([("DP", "34"), ("JQ_VS_HC_SOM", "1"), ("JQ_MT_HC_SOM", "1"), ("JQ_FOO_HC_SOM", "1"), ("JQ_VS_AF", "0.5"), ("JQ_MT_AF", "0.43"), ("JQ_MT_DP", "32"), ("JQ_VS_DP", "23"), ("{0}AF".format(consensus.JQ_CONSENSUS_TAG), "0"), ("{0}AF".format(consensus.JQ_CONSENSUS_TAG), "0.47"), ("{0}SOM_SUM".format(consensus.JQ_CONSENSUS_TAG), "3"), ("{0}DP".format(consensus.JQ_CONSENSUS_TAG), "27.5"), ('{0}AF_RANGE'.format(consensus.JQ_CONSENSUS_TAG), '0.07'), ('{0}DP_RANGE'.format(consensus.JQ_CONSENSUS_TAG), '9.0')])
        self.assertEquals(expected_dict, actual_dict)
        self.assertEquals([0.07], af)
        self.assertEquals([9.0], dp)

    def test_processLine_consensus(self):
        line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_VS_HC_SOM:DP:JQ_MT_HC_SOM:JQ_VS_AF:JQ_VS_DP\t1:0:1:0.2:43\t0:1:0:0.42:0\t1:3:0:0.0:12\n"
        new_line = process_line(line, [], "", "", [], "", "", "consensus")

        expected_line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_VS_HC_SOM:DP:JQ_MT_HC_SOM:JQ_VS_AF:JQ_VS_DP:{0}AF:{0}SOM_SUM:{0}DP:{0}AF_RANGE:{0}DP_RANGE\t1:0:1:0.2:43:0.2:2:43.0:0.0:0.0\t0:1:0:0.42:0:0.42:0:0.0:0.0:0.0\t1:3:0:0.0:12:0.0:1:12.0:0.0:0.0\n".format(consensus.JQ_CONSENSUS_TAG)
        self.assertEquals(expected_line, new_line)

    def test_processLine_consensusNoJQTags(self):
        line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tFOO:DP\t1:0\t0:2\t1:3\n"
        new_line = process_line(line, [], "", "", [], "", "", "consensus")
        expected_line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tFOO:DP:{0}AF:{0}SOM_SUM:{0}DP:{0}AF_RANGE:{0}DP_RANGE\t1:0:0.0:0:0.0:0.0:0.0\t0:2:0.0:0:0.0:0.0:0.0\t1:3:0.0:0:0.0:0.0:0.0\n".format(consensus.JQ_CONSENSUS_TAG)
        self.assertEquals(expected_line, new_line)

    def test_processLine_zscore(self):
        af_range = [0.1, 1.2, 1.0]
        af_mean = sum(af_range)/len(af_range) if af_range != [] else ""
        af_std = numpy.std(af_range) if af_range != [] else ""

        dp_range = [10.0, 1.0, 3.0]
        dp_mean = sum(dp_range)/len(dp_range) if dp_range != [] else ""
        dp_std = numpy.std(dp_range) if dp_range != [] else ""

        line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:JQ_DP_VS:JQ_AF_VS:{0}AF_RANGE:{0}DP_RANGE\t1:10.0:0.1:0.1:10.0\t0:1.0:1.2:1.2:1.0\t1:3.0:1.0:1.0:3.0\n".format(consensus.JQ_CONSENSUS_TAG)
        actual_line = process_line(line, af_range, af_mean, af_std, dp_range, dp_mean, dp_std, "zscore")

        expected_line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:JQ_DP_VS:JQ_AF_VS:{0}AF_RANGE:{0}DP_RANGE:{0}AF_RANGE_ZSCORE:{0}DP_RANGE_ZSCORE\t1:10.0:0.1:0.1:10.0:-1.39:1.38\t0:1.0:1.2:1.2:1.0:0.91:-0.95\t1:3.0:1.0:1.0:3.0:0.49:-0.43\n".format(consensus.JQ_CONSENSUS_TAG)
        self.assertEquals(expected_line, actual_line)

    def test_createConsensusDict_som(self):
        input_dict = OrderedDict([("JQ_SOM_VS", "0,1"), ("DP", "23")])
        consensus_dict = create_consensus_dict("JQ_SOM_VS", "0,1", input_dict, {}, "int")
        expected_dict = {0: 0, 1: 1}
        self.assertEquals(expected_dict, consensus_dict)

    def test_createConsensusDict_af(self):
        input_dict = OrderedDict([("JQ_AF_VS", "0.02,0.23"), ("DP", "23")])
        consensus_dict = create_consensus_dict("JQ_AF_VS", "0.02,0.23", input_dict, {}, "float")
        expected_dict = {0: 0.02, 1: 0.23}
        self.assertEquals(expected_dict, consensus_dict)

    def test_createConsensusDict_dp(self):
        input_dict = OrderedDict([("JQ_DP_VS", "2,23"), ("DP", "23")])
        consensus_dict = create_consensus_dict("JQ_DP_VS", "2,23", input_dict, {}, "float")
        expected_dict = {0: 2.0, 1: 23.0}
        self.assertEquals(expected_dict, consensus_dict)

    def test_getConsensusSom(self):
        field_dict = {0:1, 1:4, 2:0}
        consensus = get_consensus_som(field_dict)
        self.assertEquals("1,4,0", consensus)

        field_dict = {0:1}
        consensus = get_consensus_som(field_dict)
        self.assertEquals("1", consensus)

    def test_getConsensus(self):
        tags = ["1", "2"]
        field_dict = {0:1.0, 1:4.4, 2:0.1}
        consensus = get_consensus(tags, field_dict)
        self.assertEquals("0.5,2.2,0.05", consensus)

        tags = ["1"]
        field_dict = {0:1.0}
        consensus = get_consensus(tags, field_dict)
        self.assertEquals("1.0", consensus)

    def test_addZscore(self):
        meta_headers = []
        header = "#CHROM"
        lines = ["1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:JQ_DP_VS:JQ_AF_VS:{0}AF_RANGE:{0}DP_RANGE:{0}AF_RANGE_ZSCORE:{0}DP_RANGE_ZSCORE\t1:10.0:0.1:0.1:10.0:-1.39:1.38\t0:1.0:1.2:1.2:1.0:0.91:-0.95\t1:3.0:1.0:1.0:3.0:0.49:-0.43\n".format(consensus.JQ_CONSENSUS_TAG)]
        writer = MockWriter()
        output_file = "output"
        af_range = [0.1, 1.2, 1.0]
        dp_range = [10.0, 1.0, 3.0]
        add_zscore(meta_headers, header, lines, writer, output_file, af_range, dp_range)
        self.assertEquals(['##FORMAT=<ID={0}AF_RANGE_ZSCORE,Number=A,Type=Integer,Description="Jacquard measure of consistency of allele frequencies among callers = (sample AF range - population mean AF range)/standard dev(population AF range)">'.format(consensus.JQ_CONSENSUS_TAG),
                           '##jacquard.consensus.{0}AF_RANGE_ZSCORE.mean_AF_range=0.77'.format(consensus.JQ_CONSENSUS_TAG),
                            '##jacquard.consensus.{0}AF_RANGE_ZSCORE.standard_deviation=0.48'.format(consensus.JQ_CONSENSUS_TAG), 
                            '##FORMAT=<ID={0}DP_RANGE_ZSCORE,Number=A,Type=Integer,Description="Jacquard measure of consistency of depth among callers = (sample DP range - population mean DP range)/standard dev(population DP range)">'.format(consensus.JQ_CONSENSUS_TAG), 
                            '##jacquard.consensus.{0}DP_RANGE_ZSCORE.mean_DP_range=4.67'.format(consensus.JQ_CONSENSUS_TAG), 
                            '##jacquard.consensus.{0}DP_RANGE_ZSCORE.standard deviation_DP_range=3.86'.format(consensus.JQ_CONSENSUS_TAG), 
                            '#CHROM', 
                            '1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:JQ_DP_VS:JQ_AF_VS:{0}AF_RANGE:{0}DP_RANGE:{0}AF_RANGE_ZSCORE:{0}DP_RANGE_ZSCORE\t1:10.0:0.1:0.1:10.0:-1.39:1.38\t0:1.0:1.2:1.2:1.0:0.91:-0.95\t1:3.0:1.0:1.0:3.0:0.49:-0.43'.format(consensus.JQ_CONSENSUS_TAG)], writer.lines())

    def test_calculateZscore(self):
        af_range = [0.1, 1.2, 1.0]
        af_mean = sum(af_range)/len(af_range) if af_range != [] else ""
        af_std = numpy.std(af_range) if af_range != [] else ""

        dp_range = [10.0, 1.0, 3.0]
        dp_mean = sum(dp_range)/len(dp_range) if dp_range != [] else ""
        dp_std = numpy.std(dp_range) if dp_range != [] else ""

        combined_dict = {"{0}AF_RANGE".format(consensus.JQ_CONSENSUS_TAG) : "0.1", "{0}DP_RANGE".format(consensus.JQ_CONSENSUS_TAG) : "10.0"}
        combined_dict = calculate_zscore(af_mean, af_std, dp_mean, dp_std, combined_dict)

        self.assertEquals({'{0}DP_RANGE_ZSCORE'.format(consensus.JQ_CONSENSUS_TAG): '1.38', '{0}AF_RANGE_ZSCORE'.format(consensus.JQ_CONSENSUS_TAG): '-1.39', '{0}AF_RANGE'.format(consensus.JQ_CONSENSUS_TAG): '0.1', '{0}DP_RANGE'.format(consensus.JQ_CONSENSUS_TAG): '10.0'}, combined_dict)
        
    def test_functional_consensus(self):
        with TempDirectory() as output_dir:
            module_testdir = os.path.dirname(os.path.realpath(__file__))+"/functional_tests/05_consensus"
            input_dir = os.path.join(module_testdir,"input")
            args = Namespace(input=os.path.join(input_dir,os.listdir(input_dir)[0]), 
                         output=os.path.join(output_dir.path,"consensus.vcf"))
            
            execution_context = ["##jacquard.version={0}".format(utils.__version__),
                "##jacquard.command=",
                "##jacquard.cwd="]
            
            consensus.execute(args,execution_context)
            
            output_file = glob.glob(os.path.join(output_dir.path, "*.vcf"))[0]
            
            actual_file = FileReader(output_file)
            actual_file.open()
            actual = []
            for line in actual_file.read_lines():
                actual.append(line)
            actual_file.close()
            
            module_outdir = os.path.join(module_testdir,"benchmark")
            print os.listdir(module_outdir)
            output_file = os.listdir(module_outdir)[0]
            expected_file = FileReader(os.path.join(module_outdir,output_file))
            expected_file.open()
            expected = []
            for line in expected_file.read_lines():
                expected.append(line)
            expected_file.close()
            
            self.assertEquals(len(expected), len(actual))
            
            self.assertEquals(30, len(actual))
            
            for i in xrange(len(expected)):
                if expected[i].startswith("##jacquard.cwd="):
                    self.assertTrue(actual[i].startswith("##jacquard.cwd="))
                elif expected[i].startswith("##jacquard.command="):
                    self.assertTrue(actual[i].startswith("##jacquard.command="))
                else:
                    self.assertEquals(expected[i].rstrip(), actual[i].rstrip()) 

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

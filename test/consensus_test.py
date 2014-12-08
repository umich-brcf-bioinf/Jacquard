# pylint: disable=C0111,C0301,R0904,C0111,W0212
from argparse import Namespace
import glob
import os
from StringIO import StringIO
import sys
from testfixtures import TempDirectory
import unittest

import jacquard.utils as utils
import jacquard.consensus as consensus
from jacquard.variant_callers import consensus_helper
import jacquard.logger as logger
import jacquard.vcf as vcf

class MockFileWriter(object):
    def __init__(self):
        self._content = []
        self.opened = False
        self.closed = False

    def open (self):
        self.opened = True

    def write(self, content):
        if content == None:
            return
        self._content.extend(content.splitlines())

    def lines(self):
        return self._content

    def close(self):
        self.closed = True

class MockVcfReader(object):
    def __init__(self, input_filepath="vcfName", metaheaders=["##metaheaders"],
                 column_header="#header"):
        self.input_filepath = input_filepath
        self.metaheaders = metaheaders
        self.column_header = column_header
        self.opened = False
        self.closed = False

    def open(self):
        self.opened = True

    def vcf_records(self):
        return iter(["foo"])

    def close(self):
        self.closed = True

class MockConsensusTag(object):
    def __init__(self, metaheader = "##consensus_metaheader", all_ranges = []):
        self.metaheader = metaheader
        self.format_called = False
        self.all_ranges = all_ranges
        self.name = "af"

    def format(self):
        self.format_called = True

    def _roundTwoDigits(self, value):
        split_value = value.split(".")

        if len(split_value[1]) <= 2:
            if split_value[1] == '0':
                return split_value[0]
            return value

        else:
            return str(round(100 * float(value))/100)

class MockConsensusHelper(object):
    def __init__(self, tag, ranges=[0.0]):
        self.tags = [tag]
        self.add_tags_called = False
        self.add_zscore_called = False
        self.ranges = ranges

    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.format_called = True
        self.add_tags_called = True

    def add_zscore(self, vcf_record, pop_values):
        for tag in self.tags:
            tag.format_called = True
        self.add_zscore_called = True

    def get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

mock_log_called = False

def mock_log(msg, *args):
    global mock_log_called
    mock_log_called = True

class Consensus2TestCase(unittest.TestCase):
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

    def test_write_metaheaders(self):
        file_writer = MockFileWriter()
        vcf_reader = MockVcfReader()
        cons_helper = MockConsensusHelper(MockConsensusTag())
        consensus._write_metaheaders(cons_helper,
                           vcf_reader,
                           file_writer,
                           ["execution_context"])
        expected = ["##metaheaders", "execution_context",
                    "##consensus_metaheader", "#header"]
        self.assertEquals(expected,file_writer._content)

    def test_write_execution_metaheaders(self):
        pop_values = {"AF": [0.3, 0.1], "DP": [5, 18]}
        actual = consensus._get_execution_metaheaders(pop_values)
        expected = ["##jacquard.consensus.JQ_CONS_DP_RANGE.mean_DP_range=5",
                    "##jacquard.consensus.JQ_CONS_DP_ZSCORE.std_DP_range=18",
                    "##jacquard.consensus.JQ_CONS_AF_RANGE.mean_AF_range=0.3",
                    "##jacquard.consensus.JQ_CONS_AF_ZSCORE.std_AF_range=0.1"]

        self.assertEquals("\n".join(expected), actual)

    def test_add_consensus_tags(self):
        file_writer = MockFileWriter()
        vcf_reader = MockVcfReader()
        tag = MockConsensusTag()
        cons_helper = MockConsensusHelper(tag)

        consensus._add_consensus_tags(cons_helper, vcf_reader, file_writer)

        self.assertTrue(cons_helper.add_tags_called)
        self.assertTrue(tag.format_called)

    def test_add_zscore(self):
        file_writer = MockFileWriter()
        vcf_reader = MockVcfReader()
        tag = MockConsensusTag(all_ranges=[])
        ranges = [0.2, 0.3, 0.4]
        cons_helper = MockConsensusHelper(tag, ranges)
        pop_values = {"allele_freq_tag": [0.3, 0.1]}

        consensus._add_zscore(cons_helper, vcf_reader, file_writer, pop_values)

        self.assertTrue(cons_helper.add_zscore_called)
        self.assertTrue(tag.format_called)

    def test_add_zscore_rangesAreZero(self):
        file_writer = MockFileWriter()
        vcf_reader = MockVcfReader()
        tag = MockConsensusTag(all_ranges=[])
        ranges = [0.0, 0.0, 0.0]
        cons_helper = MockConsensusHelper(tag, ranges)
        pop_values = {"allele_freq_tag": [0.0, 0.0]}

        consensus._add_zscore(cons_helper, vcf_reader, file_writer, pop_values)

        self.assertTrue(cons_helper.add_zscore_called)
        self.assertTrue(tag.format_called)

    def test_add_zscore_oneCaller_NoRange(self):
        file_writer = MockFileWriter()
        vcf_reader = MockVcfReader()
        tag = MockConsensusTag(all_ranges=[])
        ranges = []
        cons_helper = MockConsensusHelper(tag, ranges)
        pop_values = {"allele_freq_tag": [0.0, 0.0]}

        consensus._add_zscore(cons_helper, vcf_reader, file_writer, pop_values)

        self.assertTrue(cons_helper.add_zscore_called)
        self.assertTrue(tag.format_called)

    def test_execute_outputFile(self):
        input_data = ("##blah\n#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|"+
                    "SAMPLE\n1|42|.|A|G|.|PASS|INFO|JQ_VS_AF:JQ_MT_AF:JQ_VS_DP:JQ_MT_DP|0.2:0.4:30:45")\
                    .replace("|","\t")
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("foo.vcf", input_data)
            input_file = os.path.join(input_dir.path,"foo.vcf")
            output_file = os.path.join(output_dir.path,"baz.vcf")
            args = Namespace(input=input_file,
                             output=output_file,
                             column_specification=None)
            consensus.execute(args,["##foo"])
            output_dir.check("baz.vcf")

    def test_execute_badInputFile(self):
        input_data = ("##blah\n#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|"+
                    "SAMPLE\n1|42|.|A|G|.|PASS|INFO|JQ_VS_AF:JQ_MT_AF:JQ_VS_DP:JQ_MT_DP|0.2:0.4:30:45")\
                    .replace("|","\t")
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("foo.txt", input_data)
            input_file = os.path.join(input_dir.path,"foo.txt")
            output_file = os.path.join(output_dir.path,"baz.vcf")
            args = Namespace(input=input_file,
                             output=output_file,
                             column_specification=None)
            self.assertRaisesRegexp(utils.JQException, "Input file "+
                                    "\[.*foo.txt\] must be a VCF file.",
                                    consensus.execute, args, ["##foo"])

    def test_execute_badOutputFile(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("foo.vcf","")
            input_file = os.path.join(input_dir.path,"foo.vcf")
            output_file = os.path.join(output_dir.path,"baz.txt")
            args = Namespace(input=input_file,
                             output=output_file,
                             column_specification=None)
            self.assertRaisesRegexp(utils.JQException, "Output file "+
                                    "\[.*baz.txt\] must be a VCF file.",
                                    consensus.execute, args, ["##foo"])

    def test_execute_outputDirectory(self):
        input_data = ("##blah\n#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|"+
                    "SAMPLE\n1|42|.|A|G|.|PASS|INFO|JQ_VS_AF:JQ_MT_AF:JQ_VS_DP:JQ_MT_DP|0.2:0.4:30:45")\
                    .replace("|","\t")
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("foo.vcf", input_data)
            input_file = os.path.join(input_dir.path,"foo.vcf")
            args = Namespace(input=input_file,
                             output=output_dir.path,
                             column_specification=None)
            consensus.execute(args,["##foo"])
            output_dir.check("consensus.vcf")

    def test_functional_consensus(self):
        file_dirname = os.path.dirname(os.path.realpath(__file__))
        module_testdir = os.path.join(file_dirname,
                                      "functional_tests",
                                      "05_consensus")
        input_dir = os.path.join(module_testdir, "input")
        output_dir = os.path.join(module_testdir, "output")

        args = Namespace(input=os.path.join(input_dir, 
                                            os.listdir(input_dir)[0]),
                         output=os.path.join(output_dir,
                                             "consensus.vcf"))

        execution_context = ["##jacquard.version={0}"\
                             .format(utils.__version__),
                             "##jacquard.command=functional test",
                             "##jacquard.cwd=foo"]

        consensus.execute(args,execution_context)

        output_file = glob.glob(os.path.join(output_dir, "*.vcf"))[0]

        actual_file = vcf.FileReader(output_file)
        actual_file.open()
        actual = []
        for line in actual_file.read_lines():
            actual.append(line)
        actual_file.close()

        module_outdir = os.path.join(module_testdir,"benchmark")

        output_file = os.listdir(module_outdir)[0]
        expected_file = vcf.FileReader(os.path.join(module_outdir,output_file))
        expected_file.open()
        expected = []
        for line in expected_file.read_lines():
            expected.append(line)
        expected_file.close()

#             self.assertEquals(len(expected), len(actual))
        self.assertEquals(34, len(actual))

        for i in xrange(len(expected)):
            if expected[i].startswith("##jacquard.cwd="):
                self.assertTrue(actual[i].startswith("##jacquard.cwd="))
            elif expected[i].startswith("##jacquard.command="):
                self.assertTrue(actual[i].startswith("##jacquard.command="))
            else:
                self.assertEquals(expected[i].rstrip(), actual[i].rstrip())

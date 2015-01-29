#pylint: disable=invalid-name, global-statement, line-too-long
#pylint: disable=too-many-public-methods, too-few-public-methods, unused-argument
from __future__ import absolute_import
from argparse import Namespace
import os
from StringIO import StringIO
import sys
from testfixtures import TempDirectory

import jacquard.utils as utils
import jacquard.consensus as consensus
import jacquard.logger as logger
import test.test_case as test_case

class MockFileWriter(object):
    def __init__(self):
        self._content = []
        self.opened = False
        self.closed = False

    def open(self):
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
    def __init__(self, input_filepath="vcfName", metaheaders=None,
                 column_header="#header"):
        self.input_filepath = input_filepath

        if metaheaders:
            self.metaheaders = metaheaders
        else:
            self.metaheaders = ["##metaheaders"]

        self.column_header = column_header
        self.opened = False
        self.closed = False

    def open(self):
        self.opened = True

    @staticmethod
    def vcf_records():
        return iter(["foo"])

    def close(self):
        self.closed = True

class MockConsensusTag(object):
    def __init__(self,
                 metaheader="##consensus_metaheader",
                 all_ranges=None):
        self.metaheader = metaheader
        self.format_called = False

        if all_ranges:
            self.all_ranges = all_ranges
        else:
            self.all_ranges = []

        self.name = "af"

    def format(self):
        self.format_called = True

    @staticmethod
    def _round_two_digits(value):
        split_value = value.split(".")

        if len(split_value[1]) <= 2:
            if split_value[1] == '0':
                return split_value[0]
            return value

        else:
            return str(round(100 * float(value))/100)

class MockConsensusHelper(object):
    def __init__(self, tag, ranges=None):
        self.tags = [tag]
        self.add_tags_called = False
        self.add_zscore_called = False

        if ranges:
            self.ranges = ranges
        else:
            self.ranges = []

    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.format_called = True
        self.add_tags_called = True

    def add_zscore(self, vcf_record, pop_values):
        for tag in self.tags:
            tag.format_called = True
        self.add_zscore_called = True

    def get_consensus_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

MOCK_LOG_CALLED = False

def mock_log(msg, *args):
    global MOCK_LOG_CALLED
    MOCK_LOG_CALLED = True

class ConsensusTestCase(test_case.JacquardBaseTestCase):
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

    @staticmethod
    def _change_mock_logger():
        global MOCK_LOG_CALLED
        MOCK_LOG_CALLED = False

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
        expected = ["##metaheaders",
                    "execution_context",
                    "##consensus_metaheader",
                    "#header"]
        self.assertEquals(expected, file_writer._content)

    def test_add_consensus_tags(self):
        file_writer = MockFileWriter()
        vcf_reader = MockVcfReader()
        tag = MockConsensusTag()
        cons_helper = MockConsensusHelper(tag)

        consensus._add_consensus_tags(cons_helper, vcf_reader, file_writer)

        self.assertTrue(cons_helper.add_tags_called)
        self.assertTrue(tag.format_called)

    def test_execute_outputFile(self):
        input_data = self.entab(\
"""##blah\n#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SAMPLE
1|42|.|A|G|.|PASS|INFO|JQ_VS_AF:JQ_MT_AF:JQ_VS_DP:JQ_MT_DP|0.2:0.4:30:45""")
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("foo.vcf", input_data)
            input_file = os.path.join(input_dir.path, "foo.vcf")
            output_file = os.path.join(output_dir.path, "baz.vcf")
            args = Namespace(input=input_file,
                             output=output_file,
                             column_specification=None)
            consensus.execute(args, ["##foo"])
            output_dir.check("baz.vcf")

    def test_execute_badInputFile(self):
        input_data = self.entab(\
"""##blah
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SAMPLE
1|42|.|A|G|.|PASS|INFO|JQ_VS_AF:JQ_MT_AF:JQ_VS_DP:JQ_MT_DP|0.2:0.4:30:45""")
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("foo.txt", input_data)
            input_file = os.path.join(input_dir.path, "foo.txt")
            output_file = os.path.join(output_dir.path, "baz.vcf")
            args = Namespace(input=input_file,
                             output=output_file,
                             column_specification=None)
            self.assertRaisesRegexp(utils.JQException,
                                    r"Input file \[.*foo.txt\] must be a VCF file.",
                                    consensus.execute, args, ["##foo"])

    def test_execute_badOutputFile(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("foo.vcf", "")
            input_file = os.path.join(input_dir.path, "foo.vcf")
            output_file = os.path.join(output_dir.path, "baz.txt")
            args = Namespace(input=input_file,
                             output=output_file,
                             column_specification=None)
            self.assertRaisesRegexp(utils.JQException, "Output file "+
                                    r"\[.*baz.txt\] must be a VCF file.",
                                    consensus.execute, args, ["##foo"])

    def test_execute_outputDirectory(self):
        input_data = self.entab(\
"""##blah
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SAMPLE
1|42|.|A|G|.|PASS|INFO|JQ_VS_AF:JQ_MT_AF:JQ_VS_DP:JQ_MT_DP|0.2:0.4:30:45""")
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("foo.vcf", input_data)
            input_file = os.path.join(input_dir.path, "foo.vcf")
            args = Namespace(input=input_file,
                             output=output_dir.path,
                             column_specification=None)
            consensus.execute(args, ["##foo"])
            output_dir.check("consensus.vcf")

class ConsensusFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_consensus(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "05_consensus")
            input_dir = os.path.join(module_testdir, "input", "tiny_strelka.merged.vcf")
            output_file = os.path.join(output_dir.path, "consensus.vcf")

            command = ["consensus", input_dir, output_file, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

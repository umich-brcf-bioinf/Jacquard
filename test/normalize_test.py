#pylint: disable=line-too-long, global-statement, unused-argument
#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
from __future__ import absolute_import

from StringIO import StringIO
from argparse import Namespace
import os
import sys
import unittest

from testfixtures import TempDirectory

import jacquard.logger as logger
import jacquard.normalize as normalize
import jacquard.utils as utils
from jacquard.variant_callers import variant_caller_factory
from jacquard.vcf import FileReader, FileWriter
import test.test_case as test_case

MOCK_LOG_CALLED = False

def mock_log(msg, *args):
    global MOCK_LOG_CALLED
    MOCK_LOG_CALLED = True

class MockCallerFactory(object):
    def __init__(self, caller):
        self.caller = caller
        self.last_filename = None

    def get_caller(self, metaheaders, column_header, name):
        self.last_filename = name
        return self.caller

class MockCaller(object):
    def __init__(self, name="MockCaller", metaheaders=None):
        self.name = name
        if metaheaders:
            self.metaheaders = metaheaders
        else:
            self.metaheaders = ["##mockMetaheader1"]
        self.file_name_search = "snps|indels"

    @staticmethod
    def add_tags(vcfRecord):
        return vcfRecord

    @staticmethod
    def decorate_files(filenames, decorator):
        return filenames[0]+"foo"

    def get_new_metaheaders(self):
        return self.metaheaders

class MockFileReader(object):
    def __init__(self, input_filepath="/foo/mockFileReader.txt", content=None):
        self.input_filepath = input_filepath
        self.file_name = os.path.basename(input_filepath)

        if content:
            self._content = content
        else:
            self._content = []

        self.open_was_called = False
        self.close_was_called = False

    def open(self):
        self.open_was_called = True

    def read_lines(self):
        for line in self._content:
            yield line

    def close(self):
        self.close_was_called = True

class MockFileWriter(object):
    def __init__(self):
        self._content = []
        self.opened = False
        self.closed = False

    def open (self):
        self.opened = True

    def write(self, content):
        self._content.extend(content.splitlines())

    def lines(self):
        return self._content

    def close(self):
        self.closed = True

class NormalizeTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
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
        unittest.TestCase.tearDown(self)
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

    def test_predict_output_valid(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.snvs.vcf","##source=strelka\n#colHeader")
            input_dir.write("A.indels.vcf","##source=strelka\n#colHeader")
            args = Namespace(input=input_dir.path,
                             output=output_dir.path)

            existing_output_files, desired_output_files = normalize._predict_output(args)
            expected_existing_output_files = []
            expected_desired_output_files = set(["A.normalized.vcf"])

            self.assertEquals(expected_existing_output_files, existing_output_files)
            self.assertEquals(expected_desired_output_files, desired_output_files)

    def test_predict_output_invalid(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.snvs.vcf","##source=strelka\n#colHeader")
            input_dir.write("A.indels.vcf","##source=strelka\n#colHeader")
            output_dir.write("A.normalized.vcf","##source=strelka\n#colHeader")
            args = Namespace(input=input_dir.path,
                             output=output_dir.path)

            existing_output_files, desired_output_files = normalize._predict_output(args)
            expected_existing_output_files = ["A.normalized.vcf"]
            expected_desired_output_files = set(["A.normalized.vcf"])

            self.assertEquals(expected_existing_output_files, existing_output_files)
            self.assertEquals(expected_desired_output_files, desired_output_files)

    def test_validate_single_caller(self):
        with TempDirectory() as input_dir:
            strelka_file1 = input_dir.write("strelka1.vcf",
                                            "##source=strelka\n#colHeader")
            strelka_file2 = input_dir.write("strelka2.vcf",
                                            "##source=strelka\n#colHeader")
            caller = normalize._validate_single_caller([strelka_file1, strelka_file2],
                                                       variant_caller_factory.get_caller)

            self.assertEquals("Strelka", caller.name)

    def test_validate_single_caller_raisesWhenDistinctCallers(self):
        with TempDirectory() as input_dir:
            strelka_file = input_dir.write("strelka1.vcf",
                                           "##source=strelka\n#colHeader")
            varscan_file = input_dir.write("varscan.vcf",
                                           "##source=VarScan2\n" +\
                                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR")
            self.assertRaises(utils.JQException,
                              normalize._validate_single_caller,
                              [strelka_file, varscan_file],
                              variant_caller_factory.get_caller)

    def test_validate_single_caller_raisesWhenUnrecognizedCaller(self):
        with TempDirectory() as input_dir:
            strelka_file = input_dir.write("strelka1.vcf",
                                           "##source=strelka\n#colHeader")
            unrecognized_file = input_dir.write("foo.vcf",
                                                "##source=foo\n#colHeader")
            self.assertRaises(utils.JQException,
                              normalize._validate_single_caller,
                              [strelka_file, unrecognized_file],
                              variant_caller_factory.get_caller)

    def test_partition_input_files(self):
        patient_to_files = {"A": ["baz.", "bar."], "B": ["blah."]}
        caller = MockCaller()
        output_dir_path = ""
        writer_to_readers = normalize._partition_input_files(patient_to_files, output_dir_path, caller)
        self.maxDiff=None
        writerA = FileWriter("baz.foo")
        readersA = [FileReader("baz."),
                    FileReader("bar.")]
        writerB = FileWriter("blah.foo")
        readersB = [FileReader("blah.")]

        self.assertEquals({writerA: readersA, writerB: readersB},
                          writer_to_readers)

    def test_determine_caller_per_directory(self):
        with TempDirectory() as input_dir:
            A = input_dir.write("A.vcf","##source=strelka\n#colHeader")
            B = input_dir.write("B.vcf","##source=strelka\n#colHeader")
            input_files = [A, B]

            mock_caller = MockCaller()
            mock_caller_factory = MockCallerFactory(mock_caller)

            caller = normalize._determine_caller_per_directory(input_files, mock_caller_factory.get_caller)

            self.assertEquals(mock_caller,caller)
            self.assertEquals(mock_caller_factory.last_filename, "A.vcf")

    def test_execute(self):
        vcf_content1 = ('''##source=strelka
##file1
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
''').replace('|', "\t")
        vcf_content2 = ('''##source=strelka
##file2
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|10|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|10|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
''').replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("P1.strelka.snvs.vcf", vcf_content1)
            input_dir.write("P1.strelka.indels.vcf", vcf_content2)
            args = Namespace(input=input_dir.path,
                             output=output_dir.path)

            normalize.execute(args, ["extra_header1", "extra_header2"])

            output_dir.check("P1.strelka.normalized.vcf")
            with open(os.path.join(output_dir.path, "P1.strelka.normalized.vcf")) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()

        self.assertEquals(8, len(actual_output_lines), "normalize output wrong number of lines")

class NormalizeFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_normalize(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "01_normalize")
            input_dir = os.path.join(module_testdir, "input")

            command = ["normalize", input_dir, output_dir.path, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

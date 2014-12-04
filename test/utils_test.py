# pylint: disable=C0103,C0301,R0903,R0904
from collections import OrderedDict
import os
from StringIO import StringIO
import subprocess
import sys
from testfixtures import TempDirectory
import unittest

from jacquard.utils import validate_directories, write_output, sort_headers, sort_data, change_pos_to_int, combine_format_values
import jacquard.utils as utils
import jacquard.logger as logger

# logger.initialize_logger("utils")
mock_log_called = False
mock_message = ""

def mock_log(msg, *args):
    global mock_log_called
    mock_log_called = True
#     print msg.format(*[str(i) for i in args])
    
class ValidateDirectoriesTestCase(unittest.TestCase):
    def setUp(self):
        self.original_info = logger.info
        self.original_error = logger.error
        self.original_warning = logger.warning
        self.original_debug = logger.debug
        self._change_mock_logger()
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output

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
        
    def test_validateDirectories_inputDirectoryDoesntExist(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/reference_files/tag_varscan_test/foo"
        output_dir = script_dir + "/reference_files/tag_varscan_test/output"

        with self.assertRaises(SystemExit) as cm:
            validate_directories(input_dir, output_dir)
        self.assertEqual(cm.exception.code, 1)
        
        global mock_log_called
        self.assertTrue(mock_log_called)
#         self.assertRegexpMatches(self.output.getvalue(),
#                                  r"Specified input directory \[.*\] does not exist.")

    def test_validateDirectories_outputDirectoryNotCreated(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.txt",
                            "##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            unwriteable_dir = os.path.join(output_dir.path, "unwriteable")
            desired_dir = os.path.join(unwriteable_dir, "bar")

            try:
                make_unwritable_dir(unwriteable_dir)
                with self.assertRaises(SystemExit) as cm:
                    validate_directories(input_dir.path, desired_dir)

            finally:
                cleanup_unwriteable_dir(unwriteable_dir)

            self.assertEqual(cm.exception.code, 1)
            global mock_log_called
            self.assertTrue(mock_log_called)
#             self.assertRegexpMatches(self.output.getvalue(),
#                                      r"Output directory \[.*\] could not be created. Check parameters and try again")

class WriteOutputTestCase(unittest.TestCase):
    def test_writeOutput(self):
        mock_writer = MockWriter()
        headers = ["#foo", "#bar"]
        actual_sorted_variants = ["123", "456"]

        write_output(mock_writer, headers, actual_sorted_variants)
        actualLines = mock_writer.lines()

        self.assertEqual("#foo", actualLines[0])
        self.assertEqual("#bar", actualLines[1])
        self.assertEqual("123", actualLines[2])
        self.assertEqual("456", actualLines[3])

class CombineFormatValuesTestCase(unittest.TestCase):
    def test_combineFormatValues(self):
        format_tags = "DP:AF:FOO"
        sample = "23:0.32:1"
        actual_dict = combine_format_values(format_tags, sample)
        expected_dict = OrderedDict([("DP", "23"), ("AF", "0.32"), ("FOO", "1")])
        self.assertEquals(expected_dict, actual_dict)


class SortTestCase(unittest.TestCase):
    def test_sort_sortHeaders(self):
        headers = ["##foo", "##bar", "#CHROM", "##baz"]
        sorted_headers = sort_headers(headers)
        expected_sorted_headers = ["##foo", "##bar", "##baz", "#CHROM"]
        self.assertEqual(expected_sorted_headers, sorted_headers)

    def test_sort_changePosToInt(self):
        split_line = ["1", "2352", "A", "G", "foo", "DP", "234"]
        line = change_pos_to_int(split_line)
        self.assertEqual([1, 2352, "A", "G", "foo", "DP", 234], line)

    def test_sort_sortData(self):
        actual_sorted_variants = ["chr1\t2352\tA\tG\tfoo\tDP\t234",
                                  "chr1\t235234\tA\tG\tfoo\tDP\t234",
                                  "chr2\t2352\tA\tG\tfoo\tDP\t234",
                                  "chr1\t2700\tA\tG\tfoo\tDP\t345",
                                  "chr10\t2352\tA\tG\tfoo\tDP\t234",
                                  "chr1\t2\tA\tG\tfoo\tDP\t234"]
        sorted_variants = sort_data(actual_sorted_variants)
        expected_sorted_variants = ["chr1\t2\tA\tG\tfoo\tDP\t234",
                                    "chr1\t2352\tA\tG\tfoo\tDP\t234",
                                    "chr1\t2700\tA\tG\tfoo\tDP\t345",
                                    "chr1\t235234\tA\tG\tfoo\tDP\t234",
                                    "chr2\t2352\tA\tG\tfoo\tDP\t234",
                                    "chr10\t2352\tA\tG\tfoo\tDP\t234"]

        self.assertEqual(expected_sorted_variants, sorted_variants)

    def test_sort_sortData_noCHR(self):
        input_variants = ["1\t2352\tA\tG\tfoo\tDP\t234",
                          "1\t235234\tA\tG\tfoo\tDP\t234",
                          "2\t2352\tA\tG\tfoo\tDP\t234",
                          "1\t2700\tA\tG\tfoo\tDP\t345",
                          "10\t2352\tA\tG\tfoo\tDP\t234",
                          "1\t2\tA\tG\tfoo\tDP\t234"]
        actual_sorted_variants = sort_data(input_variants)
        expected_variants = ["chr1\t2\tA\tG\tfoo\tDP\t234",
                             "chr1\t2352\tA\tG\tfoo\tDP\t234",
                             "chr1\t2700\tA\tG\tfoo\tDP\t345",
                             "chr1\t235234\tA\tG\tfoo\tDP\t234",
                             "chr2\t2352\tA\tG\tfoo\tDP\t234",
                             "chr10\t2352\tA\tG\tfoo\tDP\t234"]
        self.assertEqual(expected_variants, actual_sorted_variants)

    def test_sort_sortDataSamePos(self):
        all_variants = ["chr1\t2352\tA\tG\tfoo\tDP\t234", "1\t2352\tA\tGT\tfoo\tDP\t234"]
        actual_sorted_variants = sort_data(all_variants)

        expected_variants = ['chr1\t2352\tA\tG\tfoo\tDP\t234', 'chr1\t2352\tA\tGT\tfoo\tDP\t234']
        self.assertEqual(expected_variants, actual_sorted_variants)


def is_windows_os():
    return sys.platform.lower().startswith("win")

def make_unwritable_dir(unwriteable_dir):
    os.mkdir(unwriteable_dir, 0555)

    if is_windows_os():
        FNULL = open(os.devnull, 'w')
        subprocess.call("icacls {0} /deny Everyone:W".format(unwriteable_dir), stdout=FNULL, stderr=subprocess.STDOUT)

def cleanup_unwriteable_dir(unwriteable_dir):
    if is_windows_os():
        FNULL = open(os.devnull, 'w')
        subprocess.call("icacls {0} /reset /t /c".format(unwriteable_dir), stdout=FNULL, stderr=subprocess.STDOUT)
    os.rmdir(unwriteable_dir)

class MockWriter(object):
    def __init__(self):
        self._content = []
        self.wasClosed = False

    def write(self, content):
        self._content.extend(content.splitlines())

    def lines(self):
        return self._content

    def close(self):
        self.wasClosed = True


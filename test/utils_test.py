# pylint: disable=C0103,C0301,R0903,R0904,C0111
from collections import OrderedDict
import os
from StringIO import StringIO
import subprocess
import sys
from testfixtures import TempDirectory
import unittest

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
        input_dir = script_dir + "/functional_tests/utils_test/tag_varscan_test/foo"
        output_dir = script_dir + "/functional_tests/utils_test/tag_varscan_test/output"

        with self.assertRaises(SystemExit) as cm:
            utils.validate_directories(input_dir, output_dir)
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
                    utils.validate_directories(input_dir.path, desired_dir)

            finally:
                cleanup_unwriteable_dir(unwriteable_dir)

            self.assertEqual(cm.exception.code, 1)
            global mock_log_called
            self.assertTrue(mock_log_called)

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


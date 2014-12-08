# pylint: disable=C0103,C0301,R0903,R0904,W0603,W0613,C0111
from argparse import Namespace
import glob
import os
import shutil
from StringIO import StringIO
import sys
from testfixtures import TempDirectory
import unittest

import jacquard.jacquard as jacquard
import jacquard.logger as logger
import jacquard.utils as utils
import jacquard.vcf as vcf

import jacquard.tag as tag
import jacquard.normalize as normalize
import jacquard.filter_hc_somatic as filter_hc_somatic
import jacquard.merge as merge
import jacquard.consensus as consensus
import jacquard.expand as expand

import test.mock_module as mock_module

mock_create_tmp_called = False
mock_move_tmp_contents_called = False

# pylint: disable=W0603,W0613
def mock_create_temp_directory(output_dir, force=0):
    global mock_create_tmp_called
    mock_create_tmp_called = True

    if len(os.listdir(output_dir)) != 0:
        if not force:
            sys.exit(1)

# pylint: disable=W0603,W0613
def mock_move_tmp_contents_to_original(tmp_output, output_dir):
    global mock_move_tmp_contents_called
    mock_move_tmp_contents_called = True

def _change_mock_methods():
#    global mock_create_temp_directory
    jacquard._create_temp_directory = mock_create_temp_directory

#    global mock_move_tmp_contents_to_original
    jacquard._move_tmp_contents_to_original = mock_move_tmp_contents_to_original

mock_log_called = False

def mock_log(msg, *args):
    global mock_log_called
    mock_log_called = True

class JacquardTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr
        unittest.TestCase.tearDown(self)

    #TODO (cgates): Fix
    def test_gracefulErrorMessageWhenUnanticipatedProblem(self):
        with TempDirectory() as output_dir:
            mock_module.my_exception_string = "I'm feeling angry"

            with self.assertRaises(SystemExit) as exit_code:
                jacquard.dispatch([mock_module], ["mock_module", output_dir.path])

            self.assertEqual(1, exit_code.exception.code)

            actual_messages = self.output.getvalue().rstrip().split("\n")

            self.assertEquals(4, len(actual_messages))
            self.assertRegexpMatches(actual_messages[0], "Jacquard begins")
            self.assertRegexpMatches(actual_messages[1], "Saving log to")
            self.assertRegexpMatches(actual_messages[2], "I'm feeling angry")
            self.assertRegexpMatches(actual_messages[3], "Jacquard encountered an unanticipated problem.")

    def test_create_temp_directory(self):
        with TempDirectory() as test_dir:
            output_dir = test_dir.makedir("output")
            actual_tmp_dir = jacquard._create_temp_directory(output_dir)

            self.assertTrue(os.path.exists(actual_tmp_dir), "temp dir created")
            self.assertEquals(os.path.join(output_dir, "jacquard_tmp"),
                              actual_tmp_dir)


    def test_move_tmp_contents_to_original(self):
        with TempDirectory() as output_dir:
            tmp_dir = output_dir.makedir("tmp")
            with open(os.path.join(tmp_dir, "A.txt"), "w") as file_a, \
                    open(os.path.join(tmp_dir, "B.txt"), "w") as file_b:
                file_a.write("A")
                file_b.write("B")

            jacquard._move_tmp_contents_to_original(tmp_dir, output_dir.path)
            actual_files = os.listdir(output_dir.path)
            self.assertEquals(2, len(actual_files))
            self.assertEquals(["A.txt", "B.txt"], sorted(actual_files))

class JacquardTestCase_dispatchOnly(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output
        self.original_create_temp_directory = jacquard._create_temp_directory
        self.original_move_tmp_contents = jacquard._move_tmp_contents_to_original
        _change_mock_methods()

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr
        unittest.TestCase.tearDown(self)
        jacquard._create_temp_directory = self.original_create_temp_directory
        jacquard._move_tmp_contents_to_original = self.original_move_tmp_contents

    def test_dispatch(self):
        with TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            jacquard.dispatch([mock_module], ["mock_module", output_dir.path])
            self.assertTrue(mock_module.execute_called)

            self.assertTrue(mock_create_tmp_called)
            self.assertTrue(mock_move_tmp_contents_called)

    def test_dispatch_nonEmptyOutputDir(self):
        with TempDirectory() as output_dir:
            output_dir.write("file1.vcf", "foo")
            mock_module.my_exception_string = ""

            with self.assertRaises(SystemExit) as exit_code:
                jacquard.dispatch([mock_module], ["mock_module", output_dir.path])

            self.assertEqual(1, exit_code.exception.code)

    def test_dispatch_forceNonEmptyOutputDir(self):
        with TempDirectory() as output_dir:
            output_dir.write("file1.vcf", "foo")
            mock_module.my_exception_string = ""

            jacquard.dispatch([mock_module], ["mock_module", output_dir.path, "--force"])
            self.assertTrue(1 == 1, "Force does not result in premature exit.")

    def test_dispatch_done(self):
        with TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            jacquard.dispatch([mock_module], ["mock_module", output_dir.path])
            actual_messages = self.output.getvalue().rstrip().split("\n")
            self.assertRegexpMatches(actual_messages[3], "Done")

    def test_dispatch_doneWithWarnings(self):
        with TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            logger.SHOW_WARNING = True
            jacquard.dispatch([mock_module], ["mock_module", output_dir.path])
            actual_messages = self.output.getvalue().rstrip().split("\n")
            self.assertRegexpMatches(actual_messages[3], r"Done. \(See warnings above\)")
            logger.SHOW_WARNING = False

class JacquardTestCase_FunctionalTest(unittest.TestCase):
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

#TODO: fix
    def test_functional_jacquard(self):
        with TempDirectory() as output_dir:
            file_dirname = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(file_dirname,
                                          "functional_tests",
                                          "jacquard_test")

            initial_input = os.path.join(module_testdir, "input")
            normalize_output = output_dir.path
            tag_output = os.path.join(output_dir.path, "tag")
            filter_output = os.path.join(output_dir.path, "filter_hc_somatic")
            merge_output = os.path.join(output_dir.path,
                                        "merge",
                                        "merged.vcf")
            consensus_output = os.path.join(output_dir.path,
                                            "consensus",
                                            "consensus.vcf")
            expanded_output = os.path.join(output_dir.path,
                                           "expand",
                                           "expanded.vcf")

            commands = [["normalize", initial_input, normalize_output, "--force"],
                        ["tag", normalize_output, tag_output, "--force"],
                        ["filter_hc_somatic", tag_output, filter_output, "--force"],
                        ["merge", filter_output, merge_output, "--force"],
                        ["consensus", merge_output, consensus_output, "--force"],
                        ["expand", consensus_output, expanded_output, "--force"]]

            for command in commands:
                jacquard.dispatch(jacquard._SUBCOMMANDS, command)

                tool = command[0]
                output = command[2]

                if os.path.isfile(output):
                    output_file = output
                elif os.path.isdir(output):
                    single_file = os.listdir(output)[0]
                    output_file = os.path.join(output, single_file)

                actual_file = vcf.FileReader(output_file)
                actual_file.open()
                actual = []
                for line in actual_file.read_lines():
                    actual.append(line)
                actual_file.close()

                module_outdir = os.path.join(module_testdir, tool, "benchmark")
                output_file = os.listdir(module_outdir)[0]
                expected_file = vcf.FileReader(os.path.join(module_outdir, output_file))
                expected_file.open()
                expected = []

                for line in expected_file.read_lines():
                    expected.append(line)
                expected_file.close()

                self.assertEquals(len(expected), len(actual))

                for i in xrange(len(expected)):
                    if expected[i].startswith("##jacquard.cwd="):
                        self.assertTrue(actual[i].startswith("##jacquard.cwd="))
                    elif expected[i].startswith("##jacquard.command="):
                        self.assertTrue(actual[i].startswith("##jacquard.command="))
                    else:
                        self.assertEquals(expected[i].rstrip(), actual[i].rstrip())

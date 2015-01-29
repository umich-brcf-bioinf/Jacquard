# pylint: disable=line-too-long, global-statement, unused-argument, invalid-name, too-many-locals, too-many-public-methods
from __future__ import absolute_import

from StringIO import StringIO
from argparse import Namespace
import os
import sys
import unittest

from testfixtures import TempDirectory

import jacquard.jacquard as jacquard
import jacquard.logger as logger
import jacquard.utils as utils
import test.mock_module as mock_module
import test.test_case as test_case


MOCK_NOMINATE_TMP_CALLED = False
MOCK_MOVE_TEMP_CONTENTS_CALLED = False
MOCK_PREFLIGHT_CALLED = False

def mock_preflight(output, desired_output_files, command):
    global MOCK_PREFLIGHT_CALLED
    MOCK_PREFLIGHT_CALLED = True

def mock_nominate_temp_directory(output_dir):
# pylint: disable=W0603,W0613
    global MOCK_NOMINATE_TMP_CALLED
    MOCK_NOMINATE_TMP_CALLED = True

def mock_move_tmp_contents_to_original(tmp_output, output_dir):
# pylint: disable=W0603,W0613
    global MOCK_MOVE_TEMP_CONTENTS_CALLED
    MOCK_MOVE_TEMP_CONTENTS_CALLED = True

def _change_mock_methods():
#    global mock_create_temp_directory
    jacquard._nominate_temp_directory = mock_nominate_temp_directory

#    global mock_move_tmp_contents_to_original
    jacquard._move_tmp_contents_to_original = mock_move_tmp_contents_to_original
    
#    global mock_preflight
    jacquard._preflight = mock_preflight

MOCK_LOG_CALLED = False

def mock_log(msg, *args):
    global MOCK_LOG_CALLED
    MOCK_LOG_CALLED = True

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

    def test_gracefulErrorMessageWhenUnanticipatedProblem(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            mock_module.my_exception_string = "I'm feeling angry"

            with self.assertRaises(SystemExit) as exit_code:
                jacquard.dispatch([mock_module], ["mock_module",
                                                  input_dir.path,
                                                  output_dir.path])
            self.assertEqual(1, exit_code.exception.code)

            actual_messages = self.output.getvalue().rstrip().split("\n")

            print actual_messages
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

    def test_preflight_valid(self):
        with TempDirectory() as output_dir:
            desired_output_files = set(["foo.vcf", "bar.vcf"])
            jacquard._preflight_old(output_dir.path, desired_output_files, "foo_command")
            self.assertTrue(1==1)

    def test_preflight_invalid(self):
        with TempDirectory() as output_dir:
            output_dir.write("foo.vcf","##source=strelka\n#colHeader")
            desired_output_files = set(["foo.vcf", "bar.vcf"])

            self.assertRaisesRegexp(utils.JQException,
                                r"The command \[foo_command\] would overwrite existing files \['foo.vcf'\]",
                                jacquard._preflight_old,
                                output_dir.path,
                                desired_output_files,
                                "foo_command")

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
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            jacquard.dispatch([mock_module], ["mock_module",
                                              input_dir.path,
                                              output_dir.path])
            self.assertTrue(mock_module.execute_called)
            self.assertTrue(mock_module.report_called)

            self.assertTrue(MOCK_NOMINATE_TMP_CALLED)
            self.assertTrue(MOCK_MOVE_TEMP_CONTENTS_CALLED)

    def xtest_dispatch_nonEmptyOutputDir(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            output_dir.write("file1.vcf", "foo")
            output_file = os.path.join(output_dir.path,"file1.vcf")
            mock_module.my_exception_string = ""

            with self.assertRaises(SystemExit) as exit_code:
                jacquard.dispatch([mock_module], ["mock_module",
                                                  input_dir.path,
                                                  output_file])

            self.assertEqual(1, exit_code.exception.code)

    def test_dispatch_forceNonEmptyOutputDir(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            output_dir.write("file1.vcf", "foo")
            mock_module.my_exception_string = ""

            args = Namespace(subparser_name = "mock_module",
                             input = input_dir.path,
                             output = output_dir.path,
                             force = 1)
            jacquard.dispatch([mock_module], ["mock_module",
                                              input_dir.path,
                                              output_dir.path,
                                              "--force"])

#             jacquard.dispatch([mock_module], args)
            self.assertTrue(1 == 1, "Force does not result in premature exit.")

    def test_dispatch_done(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            jacquard.dispatch([mock_module], ["mock_module",
                                              input_dir.path,
                                              output_dir.path])
            actual_messages = self.output.getvalue().rstrip().split("\n")
            self.assertRegexpMatches(actual_messages[3], "Done")

    def test_dispatch_doneWithWarnings(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            logger.SHOW_WARNING = True
            jacquard.dispatch([mock_module], ["mock_module",
                                              input_dir.path,
                                              output_dir.path])
            actual_messages = self.output.getvalue().rstrip().split("\n")
            self.assertRegexpMatches(actual_messages[3], r"Done. \(See warnings above\)")
            logger.SHOW_WARNING = False

class JacquardFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_functional_jacquard(self):
        with TempDirectory() as output_dir:
            file_dirname = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(file_dirname,
                                          "functional_tests",
                                          "jacquard_test")

            initial_input = os.path.join(module_testdir, "input")

            vs_normalize_output = os.path.join(output_dir.path, "varscan")
            sk_normalize_output = os.path.join(output_dir.path, "strelka")
            mt_normalize_output = os.path.join(output_dir.path, "mutect")

            normalize_output = os.path.join(output_dir.path, "normalize")
            tag_output = os.path.join(output_dir.path, "tag")
            filter_output = os.path.join(output_dir.path, "filter_hc_somatic")
            merge_output = os.path.join(output_dir.path, "merge2", "merged.vcf")
            consensus_output = os.path.join(output_dir.path, "consensus", "consensus.vcf")
            expanded_output = os.path.join(output_dir.path, "expand", "expanded.tsv")

            commands = [["normalize", os.path.join(initial_input, "varscan"), vs_normalize_output, "--force"],
                        ["normalize", os.path.join(initial_input, "strelka"), sk_normalize_output, "--force"],
                        ["normalize", os.path.join(initial_input, "mutect"), mt_normalize_output, "--force"]]

            for command in commands:
                benchmark_dir = os.path.join("benchmark", os.path.basename(command[1]))
                expected_dir = os.path.join(module_testdir, command[0], benchmark_dir)
                self.assertCommand(command, expected_dir)

            self.move_files([vs_normalize_output,
                             sk_normalize_output,
                             mt_normalize_output],
                            normalize_output)

            commands = [["tag", normalize_output, tag_output, "--force"],
                        ["filter_hc_somatic", tag_output, filter_output, "--force"],
                        ["merge2", filter_output, merge_output, "--force"],
                        ["consensus", merge_output, consensus_output, "--force"],
                        ["expand", consensus_output, expanded_output, "--force"]]

            for command in commands:
                expected_dir = os.path.join(module_testdir, command[0], "benchmark")
                self.assertCommand(command, expected_dir)


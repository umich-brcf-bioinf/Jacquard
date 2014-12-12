# pylint: disable=C0103,C0301,R0903,R0904,W0603,W0613,C0111,W0212
import os
import shutil
from StringIO import StringIO
import sys
from testfixtures import TempDirectory
import unittest

import jacquard.jacquard as jacquard
import jacquard.logger as logger
import jacquard.vcf as vcf
import test_case as test_case


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

class JacquardTestCase_FunctionalTest(test_case.JacquardBaseTestCase):
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
            merge_output = os.path.join(output_dir.path,
                                        "merge",
                                        "merged.vcf")
            consensus_output = os.path.join(output_dir.path,
                                            "consensus",
                                            "consensus.vcf")
            expanded_output = os.path.join(output_dir.path,
                                           "expand",
                                           "expanded.vcf")
            commands = [["normalize", os.path.join(initial_input, "varscan"), vs_normalize_output, "--force"],
                        ["normalize", os.path.join(initial_input, "strelka"), sk_normalize_output, "--force"],
                        ["normalize", os.path.join(initial_input, "mutect"), mt_normalize_output, "--force"]]

            for command in commands:
                benchmark_dir = os.path.join("benchmark", os.path.basename(command[1]))
                self.assertCommand(command, module_testdir, benchmark_dir)

            self.move_files([vs_normalize_output,
                             sk_normalize_output,
                             mt_normalize_output],
                            normalize_output)


            commands = [["tag", normalize_output, tag_output, "--force"],
                        ["filter_hc_somatic", tag_output, filter_output, "--force"],
                        ["merge", filter_output, merge_output, "--force"],
                        ["consensus", merge_output, consensus_output, "--force"],
                        ["expand", consensus_output, expanded_output, "--force"]]

            for command in commands:
                self.assertCommand(command, module_testdir, "benchmark")


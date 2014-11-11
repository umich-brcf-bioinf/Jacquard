# pylint: disable=C0103,C0301,R0903,R0904,W0603,W0613
import os
from StringIO import StringIO
import sys
from testfixtures import TempDirectory
import unittest

import jacquard.jacquard as jacquard
import test.mock_module as mock_module

mock_create_tmp_called = False
mock_move_tmp_contents_called = False

# pylint: disable=W0603,W0613
def mock_create_temp_directory(output_dir):
    global mock_create_tmp_called
    mock_create_tmp_called = True

# pylint: disable=W0603,W0613
def mock_move_tmp_contents_to_original(tmp_output, output_dir):
    global mock_move_tmp_contents_called
    mock_move_tmp_contents_called = True

def _change_mock_methods():
#    global mock_create_temp_directory
    jacquard._create_temp_directory = mock_create_temp_directory

#    global mock_move_tmp_contents_to_original
    jacquard._move_tmp_contents_to_original = mock_move_tmp_contents_to_original


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
    def Xtest_gracefulErrorMessageWhenUnanticipatedProblem(self):
        with TempDirectory() as output_dir:
            mock_module.my_exception_string = "I'm feeling angry"

            with self.assertRaises(SystemExit) as exit_code:
                jacquard.dispatch([mock_module], ["mock_module", output_dir.path])

            self.assertEqual(1, exit_code.exception.code)

            actual_messages = self.output.getvalue().rstrip().split("\n")

            self.assertEquals(5, len(actual_messages))
            self.assertRegexpMatches(actual_messages[0], "Jacquard begins")
            self.assertRegexpMatches(actual_messages[1], "Saving log to")
            self.assertRegexpMatches(actual_messages[3], "I'm feeling angry")
            self.assertRegexpMatches(actual_messages[4], "Jacquard encountered an unanticipated problem.")

    def test_create_temp_directory(self):
        with TempDirectory() as output_dir:
            actual_tmp_dir = jacquard._create_temp_directory(output_dir.path)
            self.assertTrue(os.path.exists(actual_tmp_dir), "temp dir created")
            self.assertEquals(os.path.join(output_dir.path, "tmp"),
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
            self.assertEquals(["A.txt", "B.txt"], actual_files)

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

    #TODO (cgates): Fix
    def Xtest_dispatch(self):
        with TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            jacquard.dispatch([mock_module], ["mock_module", output_dir.path])
            self.assertTrue(mock_module.execute_called)

#            global mock_create_tmp_called
            self.assertTrue(mock_create_tmp_called)
#            global mock_move_tmp_contents_called
            self.assertTrue(mock_move_tmp_contents_called)

    #TODO (cgates): Fix
    def Xtest_dispatch_nonEmptyOutputDir(self):
        with TempDirectory() as output_dir:
            output_dir.write("file1.vcf", "foo")
            mock_module.my_exception_string = ""

            with self.assertRaises(SystemExit) as exit_code:
                jacquard.dispatch([mock_module], ["mock_module", output_dir.path])

            self.assertEqual(1, exit_code.exception.code)

    #TODO (cgates): Fix
    def test_dispatch_forceNonEmptyOutputDir(self):
        with TempDirectory() as output_dir:
            output_dir.write("file1.vcf", "foo")
            mock_module.my_exception_string = ""

            jacquard.dispatch([mock_module], ["mock_module", output_dir.path, "--force"])
            self.assertTrue(1 == 1, "Force does not result in premature exit.")


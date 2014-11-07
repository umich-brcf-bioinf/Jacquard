# pylint: disable=C0103,C0301,R0903,R0904
import os
from StringIO import StringIO
import sys
from testfixtures import TempDirectory
import unittest

import jacquard.jacquard as jacquard
import jacquard.utils as utils
import test.mock_module as mock_module

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

mock_create_tmp_called = False
mock_move_tmp_contents_called = False

def mock_create_temp_directory(output_dir):
    global mock_create_tmp_called
    mock_create_tmp_called = True
    
def mock_move_tmp_contents_to_original(tmp_output, output_dir):
    global mock_move_tmp_contents_called
    mock_move_tmp_contents_called = True

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
        output_dir = os.path.join(TEST_DIRECTORY, "reference_files", "jacquard_test", "output_dir")
        tmp_dir = jacquard._create_temp_directory(output_dir)
        self.assertEquals(os.path.join(output_dir, "tmp"), tmp_dir)

        os.rmdir(tmp_dir)

    def test_move_tmp_contents_to_original(self):
        output_dir = os.path.join(TEST_DIRECTORY, "reference_files", "jacquard_test", "output_dir")
        tmp_dir = os.path.join(output_dir, "tmp")

        try:
            os.mkdir(tmp_dir)
        except:
            pass

        open(os.path.join(tmp_dir, "foo.txt"), "a").close()
        open(os.path.join(tmp_dir, "bar.txt"), "a").close()

        jacquard._move_tmp_contents_to_original(tmp_dir, output_dir)

        self.assertEquals(2, len(os.listdir(output_dir)))
        self.assertEquals(["bar.txt", "foo.txt"], os.listdir(output_dir))

        os.remove(os.path.join(output_dir, "foo.txt"))
        os.remove(os.path.join(output_dir, "bar.txt"))

class JacquardTestCase_dispatchOnly(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output
        self.original_create_temp_directory = jacquard._create_temp_directory
        self.original_move_tmp_contents = jacquard._move_tmp_contents_to_original
        self._change_mock_methods()

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr
        unittest.TestCase.tearDown(self)
        jacquard._create_temp_directory = self.original_create_temp_directory
        jacquard._move_tmp_contents_to_original = self.original_move_tmp_contents

    def _change_mock_methods(self):
        global mock_create_temp_directory
        jacquard._create_temp_directory = mock_create_temp_directory

        global mock_move_tmp_contents_to_original
        jacquard._move_tmp_contents_to_original = mock_move_tmp_contents_to_original

    def test_dispatch(self):
        with TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            jacquard.dispatch([mock_module], ["mock_module", output_dir.path])
            self.assertTrue( mock_module.execute_called)

            global mock_create_tmp_called
            self.assertTrue(mock_create_tmp_called)
            global mock_move_tmp_contents_called
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

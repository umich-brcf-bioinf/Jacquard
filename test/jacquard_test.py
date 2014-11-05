# pylint: disable=C0103,C0301,R0903,R0904
import unittest
from StringIO import StringIO
import sys
from testfixtures import TempDirectory

import jacquard.jacquard as jacquard
import test.mock_module as mock_module

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

    def test_dispatch(self):
        with TempDirectory() as output_dir:
            jacquard.dispatch([mock_module], ["mock_module", output_dir.path])
            self.assertEqual(True, mock_module.execute_called)
        
    def test_gracefulErrorMessageWhenUnanticipatedProblem(self):
        with TempDirectory() as output_dir:
            mock_module.my_exception_string = "I'm feeling angry" 

            with self.assertRaises(SystemExit) as exit_code:
                jacquard.dispatch([mock_module], ["mock_module", output_dir.path])
    
            self.assertEqual(1, exit_code.exception.code)
    
            actual_messages = self.output.getvalue().rstrip().split("\n")
            print actual_messages
            self.assertEquals(5, len(actual_messages))
            self.assertRegexpMatches(actual_messages[0], "Jacquard begins")
            self.assertRegexpMatches(actual_messages[1], "Saving log to")
            self.assertRegexpMatches(actual_messages[3], "I'm feeling angry")
            self.assertRegexpMatches(actual_messages[4], "Jacquard encountered an unanticipated problem.")


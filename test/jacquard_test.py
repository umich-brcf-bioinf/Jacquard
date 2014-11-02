# pylint: disable=C0103,C0301,R0903,R0904
import unittest
from StringIO import StringIO
import sys

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
        jacquard.dispatch([mock_module], ["mock_module"])
        self.assertEqual(True, mock_module.execute_called)

    def test_gracefulErrorMessageWhenUnanticipatedProblem(self):
        mock_module.my_exception_string = "I'm feeling angry"

        with self.assertRaises(SystemExit) as exit_code:
            jacquard.dispatch([mock_module], ["mock_module"])

        self.assertEqual(1, exit_code.exception.code)

        actual_messages = self.output.getvalue().rstrip().split("\n")
        self.assertEquals(4, len(actual_messages))
        self.assertRegexpMatches(actual_messages[0], "Jacquard begins")
        self.assertRegexpMatches(actual_messages[1], "Saving log to")
        self.assertRegexpMatches(actual_messages[2], "I'm feeling angry")
        self.assertRegexpMatches(actual_messages[3], "Jacquard encountered an unanticipated problem.")


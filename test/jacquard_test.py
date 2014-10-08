# pylint: disable=R0904
import unittest
import jacquard.jacquard as jacquard
import test.mock_module as mock_module


class JacquardTestCase(unittest.TestCase):
    def test_main(self):
        jacquard.dispatch([mock_module], ["mock_module"])
        self.assertEqual(True, mock_module.execute_called)

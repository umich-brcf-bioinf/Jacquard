# pylint: disable=R0904
import unittest
from jacquard.jacquard import dispatch
import test.mock_module as mock_module

class JacquardTestCase(unittest.TestCase):
    def test_main(self):
        dispatch([mock_module], ["test.mock_module"])
        self.assertEqual(True, mock_module.execute_called)

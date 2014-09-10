import argparse
import unittest
from jacquard.jacquard import main
import mock_module

class JacquardTestCase(unittest.TestCase):
    def test_main(self):
        main([mock_module], ["test.mock_module"])
        self.assertEqual(True, mock_module.execute_called)
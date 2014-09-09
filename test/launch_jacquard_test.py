import subprocess
import os
import unittest

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

class LaunchJacquardTest(unittest.TestCase):
    def test_launchJacquard(self):
        p = subprocess.check_output([TEST_DIRECTORY + "/../bin/jacquard", "tag", "--help"], shell=True)
        self.assertTrue("usage: jacquard.py tag [-h] input_dir output_dir" in p)
        

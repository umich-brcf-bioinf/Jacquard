# pylint: disable=R0904,C0103
import subprocess
import os
import unittest

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

#This should be install to a tmp dir and launch
##python setup.py install --install-base=/tmp/python

class LaunchJacquardTest(unittest.TestCase):
    def Xtest_launchJacquard(self):
        process_output = subprocess.check_output(
            TEST_DIRECTORY + '/../bin/jacquard tag --help', shell=True)
        self.assertIn("usage: jacquard tag [-h] input_dir output_dir",
                      process_output)


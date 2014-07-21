#!/usr/bin/python2.7
import os
import unittest
import subprocess
import sys
import testfixtures
from testfixtures import TempDirectory
from bin.jacquard_utils import validate_directories

class ValidateDirectoriesTestCase(unittest.TestCase):
    def test_validateDirectories_inputDirectoryDoesntExist(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/tag_varscan_test/foo"
        output_dir = script_dir + "/tag_varscan_test/output"
         
        with self.assertRaises(SystemExit) as cm:
            validate_directories(input_dir, output_dir)
        self.assertEqual(cm.exception.code, 1)
     
    def test_validateDirectories_inputDirectoryUnreadable(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/tag_varscan_test/unreadable"
        output_dir = script_dir + "/tag_varscan_test/output"
 
        with self.assertRaises(SystemExit) as cm:
            validate_directories(input_dir, output_dir)
        self.assertEqual(cm.exception.code, 1)
         
    def test_validateDirectories_outputDirectoryNotCreated(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.txt","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            unwriteable_dir = os.path.join(output_dir.path,"unwriteable")
            desired_dir = os.path.join(unwriteable_dir, "bar")
 
            try:
                make_unwritable_dir(unwriteable_dir)
                 
                with self.assertRaises(SystemExit) as cm:
                    validate_directories(input_dir.path, desired_dir)

            finally:
                cleanup_unwriteable_dir(unwriteable_dir)
            
            self.assertEqual(cm.exception.code, 1)
 
def is_windows_os():
    return sys.platform.lower().startswith("win")
             
 
def make_unwritable_dir(unwriteable_dir):
    os.mkdir(unwriteable_dir, 0555)            
     
    if is_windows_os():
        FNULL = open(os.devnull, 'w')
        subprocess.call("icacls {0} /deny Everyone:W".format(unwriteable_dir), stdout=FNULL, stderr=subprocess.STDOUT)
 
def cleanup_unwriteable_dir(unwriteable_dir):
    if is_windows_os():
        FNULL = open(os.devnull, 'w')
        subprocess.call("icacls {0} /reset /t /c".format(unwriteable_dir), stdout=FNULL, stderr=subprocess.STDOUT)
    os.rmdir(unwriteable_dir)

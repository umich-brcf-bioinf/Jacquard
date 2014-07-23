#!/usr/bin/python2.7
import os
import unittest
import subprocess
import sys
import testfixtures
from testfixtures import TempDirectory
from bin.jacquard_utils import validate_directories, write_output, sort_headers, sort_data, change_pos_to_int

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
 
class WriteOutputTestCase(unittest.TestCase):
    def test_writeOutput(self):
        mock_writer = MockWriter()
        headers = ["#foo", "#bar"]
        variants =["123", "456"]
         
        write_output(mock_writer, headers, variants)
        actualLines = mock_writer.lines()
         
        self.assertEqual("#foo", actualLines[0])
        self.assertEqual("#bar", actualLines[1])
        self.assertEqual("123", actualLines[2])
        self.assertEqual("456", actualLines[3])
        
class SortTestCase(unittest.TestCase):
    def test_sort_sortHeaders(self):
        headers = ["##foo", "##bar", "#CHROM", "##baz"]
        
        sorted_headers = sort_headers(headers)
        expected_sorted_headers = ["##foo", "##bar", "##baz", "#CHROM"]
        
        self.assertEqual(expected_sorted_headers, sorted_headers)
        
    def test_sort_changePosToInt(self):
         split_line = ["1", "2352","A","G","foo","DP", "234"]
         line = change_pos_to_int(split_line)
         self.assertEqual([1,2352,"A","G","foo","DP",234], line)
         
    def test_sort_changePosToInt(self):
         split_line = ["chr1", "2352","A","G","foo","DP", "234"]
         line = change_pos_to_int(split_line)
         self.assertEqual([1,2352,"A","G","foo","DP",234], line)
     
    def test_sort_sortData(self):
        variants = ["chr1\t2352\tA\tG\tfoo\tDP\t234","chr1\t235234\tA\tG\tfoo\tDP\t234","chr2\t2352\tA\tG\tfoo\tDP\t234","chr1\t2700\tA\tG\tfoo\tDP\t345","chr10\t2352\tA\tG\tfoo\tDP\t234","chr1\t2\tA\tG\tfoo\tDP\t234"]
        
        sorted_variants = sort_data(variants)
        expected_sorted_variants = ["chr1\t2\tA\tG\tfoo\tDP\t234","chr1\t2352\tA\tG\tfoo\tDP\t234","chr1\t2700\tA\tG\tfoo\tDP\t345","chr1\t235234\tA\tG\tfoo\tDP\t234","chr2\t2352\tA\tG\tfoo\tDP\t234","chr10\t2352\tA\tG\tfoo\tDP\t234"]
        
        self.assertEqual(expected_sorted_variants, sorted_variants)
        
    def test_sort_sortData_noCHR(self):
         all_variants = ["1\t2352\tA\tG\tfoo\tDP\t234","1\t235234\tA\tG\tfoo\tDP\t234","2\t2352\tA\tG\tfoo\tDP\t234","1\t2700\tA\tG\tfoo\tDP\t345","10\t2352\tA\tG\tfoo\tDP\t234","1\t2\tA\tG\tfoo\tDP\t234"]
         variants = sort_data(all_variants)
 
         expected_variants = ["chr1\t2\tA\tG\tfoo\tDP\t234","chr1\t2352\tA\tG\tfoo\tDP\t234","chr1\t2700\tA\tG\tfoo\tDP\t345","chr1\t235234\tA\tG\tfoo\tDP\t234","chr2\t2352\tA\tG\tfoo\tDP\t234","chr10\t2352\tA\tG\tfoo\tDP\t234"]
         self.assertEqual(expected_variants, variants)
         
    def test_sort_sortDataSamePos(self):
         all_variants = ["chr1\t2352\tA\tG\tfoo\tDP\t234", "1\t2352\tA\tGT\tfoo\tDP\t234"]
         variants = sort_data(all_variants)

         expected_variants = ['chr1\t2352\tA\tG\tfoo\tDP\t234', 'chr1\t2352\tA\tGT\tfoo\tDP\t234']
         self.assertEqual(expected_variants, variants)
        
        
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
    
class MockWriter():
    def __init__(self):
        self._content = []
        self.wasClosed = False

    def write(self, content):
        self._content.extend(content.splitlines())
        
    def lines(self):
        return self._content

    def close(self):
        self.wasClosed = True


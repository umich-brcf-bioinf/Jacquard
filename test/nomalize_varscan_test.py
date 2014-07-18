#!/usr/bin/python2.7
from collections import defaultdict
import os
import unittest
import testfixtures
from testfixtures import TempDirectory
from bin.normalize_varscan import identify_merge_candidates, validate_directories, get_headers, merge_data, sort_data, change_pos_to_int, write_output, validate_split_line
class IdentifyMergeCandidatesTestCase(unittest.TestCase):
    def test_identifyMergeCandidates(self):
        in_files = ["foo.snp.vcf", "foo.indel.vcf", "bar.snp.vcf", "bar.indel.vcf", "foo.1.snp.vcf", "foo.1.indel.vcf"]
        merge_candidates = identify_merge_candidates(in_files)
        
        self.assertEqual(["foo.1.merged.vcf", "bar.merged.vcf", "foo.merged.vcf"], merge_candidates.keys())
        self.assertEqual([["foo.1.snp.vcf", "foo.1.indel.vcf"], ["bar.snp.vcf", "bar.indel.vcf"], ["foo.snp.vcf", "foo.indel.vcf"]], merge_candidates.values())

class MergeTestCase(unittest.TestCase):
    def test_merge_getHeaders(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf","##source=VarScan2\n##foobarbaz\n#CHROM\tNORMAL\tTUMOR\n123\n456\n")

            file = os.path.join(input_dir.path, "A.vcf")
            headers = get_headers(file)
            self.assertEqual(["##source=VarScan2\n", "##foobarbaz\n", "#CHROM\tNORMAL\tTUMOR\n"], headers)
            
        input_dir.cleanup()
    
    def test_merge_mergeData(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.snp.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tG\tfoo\tDP\t234\n1\t235234\tA\tG\tfoo\tDP\t234\n2\t2352\tA\tG\tfoo\tDP\t234\n")
            input_dir.write("A.indel.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2700\tA\tG\tfoo\tDP\t345\n10\t2352\tA\tG\tfoo\tDP\t234\n1\t2\tA\tG\tfoo\tDP\t234\n")

            file1 = os.path.join(input_dir.path, "A.snp.vcf")
            file2 = os.path.join(input_dir.path, "A.indel.vcf")
            all_variants = merge_data([file1, file2])

            self.assertEqual([1,2352,"A","G","foo","DP",234], all_variants[0])
            self.assertEqual([1,235234,"A","G","foo","DP",234], all_variants[1])
            self.assertEqual([2,2352,"A","G","foo","DP",234], all_variants[2])
            self.assertEqual([1,2700,"A","G","foo","DP",345], all_variants[3])
            self.assertEqual([10,2352,"A","G","foo","DP",234], all_variants[4])
            self.assertEqual([1,2,"A","G","foo","DP",234], all_variants[5])
        input_dir.cleanup()
        
    def test_merge_mergeDataSamePos(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.snp.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tG\tfoo\tDP\t234\n")
            input_dir.write("A.indel.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tGT\tfoo\tDP\t234\n")

            file1 = os.path.join(input_dir.path, "A.snp.vcf")
            file2 = os.path.join(input_dir.path, "A.indel.vcf")
            all_variants = merge_data([file1, file2])

            self.assertEqual([1,2352,"A","G","foo","DP",234], all_variants[0])
            self.assertEqual([1,2352,"A","GT","foo","DP",234], all_variants[1])
        input_dir.cleanup()
    
    def test_validateSplitLine_valid(self):
        split_line = ["chr1", "234", ".", "A", "G", ".", "PASS", "SS=5", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(0, invalid)
        
        split_line = ["chr1", "234", ".", "A", "G", ".", "PASS", "SS=2", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(0, invalid)
        
    def test_validateSplitLine_invalidAltOkaySS(self):
        split_line = ["chr1", "234", ".", "A", "G/C", ".", "PASS", "SS=5", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(0, invalid)
        
        split_line = ["chr1", "234", ".", "A", "-G", ".", "PASS", "SS=5", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(0, invalid)
        
        split_line = ["chr1", "234", ".", "A", "+G", ".", "PASS", "SS=5", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(0, invalid)
        
    def test_validateSplitLine_invalidRefOkaySS(self):
        split_line = ["chr1", "234", ".", "A/C", "G", ".", "PASS", "SS=5", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(0, invalid)
        
        split_line = ["chr1", "234", ".", "-A", "G", ".", "PASS", "SS=5", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(0, invalid)
        
        split_line = ["chr1", "234", ".", "+A", "G", ".", "PASS", "SS=5", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(0, invalid)
        
    def test_validateSplitLine_invalidAltBadSS(self):
        split_line = ["chr1", "234", ".", "A", "G/C", ".", "PASS", "SS=2", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(1, invalid)
        
        split_line = ["chr1", "234", ".", "A", "+G", ".", "PASS", "SS=2", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(1, invalid)
        
        split_line = ["chr1", "234", ".", "A", "-G", ".", "PASS", "SS=2", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(1, invalid)
        
    def test_validateSplitLine_invalidRefBadSS(self):
        split_line = ["chr1", "234", ".", "A/C", "G", ".", "PASS", "SS=2", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(1, invalid)
        
        split_line = ["chr1", "234", ".", "-A", "G", ".", "PASS", "SS=2", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(1, invalid)
        
        split_line = ["chr1", "234", ".", "+A", "G", ".", "PASS", "SS=2", "DP", "345"]
        invalid = validate_split_line(split_line, 0)
        self.assertEqual(1, invalid)
        
    def test_merge_changePosToInt(self):
        split_line = ["1", "2352","A","G","foo","DP", "234"]
        line = change_pos_to_int(split_line)
        self.assertEqual([1,2352,"A","G","foo","DP",234], line)
    
    def test_merge_sortData(self):
        all_variants = [[1,2352,"A","G","foo","DP",234],[1,235234,"A","G","foo","DP",234], [2,2352,"A","G","foo","DP",234],[1,2700,"A","G","foo","DP",345],[10,2352,"A","G","foo","DP",234],[1,2,"A","G","foo","DP",234]]
        variants = sort_data(all_variants)

        expected_variants = ["1\t2\tA\tG\tfoo\tDP\t234","1\t2352\tA\tG\tfoo\tDP\t234","1\t2700\tA\tG\tfoo\tDP\t345","1\t235234\tA\tG\tfoo\tDP\t234","2\t2352\tA\tG\tfoo\tDP\t234","10\t2352\tA\tG\tfoo\tDP\t234"]
        self.assertEqual(expected_variants, variants)
        
    def test_merge_sortDataSamePos(self):
        all_variants = [[1,2352,"A","G","foo","DP",234], [1,2352,"A","GT","foo","DP",234]]
        variants = sort_data(all_variants)

        expected_variants = ["1\t2352\tA\tG\tfoo\tDP\t234", "1\t2352\tA\tGT\tfoo\tDP\t234"]
        self.assertEqual(expected_variants, variants)
        
    def test_merge_writeOutput(self):
        mock_writer = MockWriter()
        headers = ["#foo", "#bar"]
        variants =["123", "456"]
        
        write_output(mock_writer, headers, variants)
        actualLines = mock_writer.lines()
        
        self.assertEqual("#foo", actualLines[0])
        self.assertEqual("#bar", actualLines[1])
        self.assertEqual("123", actualLines[2])
        self.assertEqual("456", actualLines[3])

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
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/tag_varscan_test/input"
        first_out_dir = script_dir + "/tag_varscan_test/unwriteable"
        
        with self.assertRaises(SystemExit) as cm:
            validate_directories(input_dir, first_out_dir + "/foo")
        self.assertEqual(cm.exception.code, 1)
        

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

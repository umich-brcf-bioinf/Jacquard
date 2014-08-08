#!/usr/bin/python2.7
from collections import defaultdict
import os
import unittest
import subprocess
import sys
import testfixtures
from testfixtures import TempDirectory
from bin.normalize_utils import VarScan, Strelka, Unknown, identify_merge_candidates, get_headers, merge_data, validate_split_line, identify_hc_variants, mark_hc_variants
        
class IdentifyMergeCandidatesTestCase(unittest.TestCase):
    def test_indentifyMergeCandidates_HC(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/normalize_strelka_test/input/"
        in_files = [input_dir + "tiny_strelka.indels.vcf", input_dir + "tiny_strelka.snvs.vcf"]
        output_dir = script_dir + "/normalize_strelka_test/output/"
        merge_candidates, hc_candidates = identify_merge_candidates(in_files, output_dir, Strelka())
        
        self.assertEqual([output_dir + "tiny_strelka.merged.vcf"], merge_candidates.keys())
        self.assertEqual([[input_dir + "tiny_strelka.indels.vcf", input_dir + "tiny_strelka.snvs.vcf"]], merge_candidates.values())
        
        
class MergeTestCase(unittest.TestCase):
    def test_merge_getHeaders(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf","##source=strelka\n##foobarbaz\n#CHROM\tNORMAL\tTUMOR\n123\n456\n")
 
            file = os.path.join(input_dir.path, "A.vcf")
            meta_headers, header = get_headers(file)
            self.assertEqual(["##source=strelka\n", "##foobarbaz\n"], meta_headers)
            self.assertEqual("#CHROM\tNORMAL\tTUMOR\n", header)
             
        input_dir.cleanup()
#     
    def test_merge_mergeData(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.snvs.vcf","##source=strelka\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tG\tfoo\tDP\t234\n1\t235234\tA\tG\tfoo\tDP\t234\n2\t2352\tA\tG\tfoo\tDP\t234\n")
            input_dir.write("A.indels.vcf","##source=strelka\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2700\tA\tG\tfoo\tDP\t345\n10\t2352\tA\tG\tfoo\tDP\t234\n1\t2\tA\tG\tfoo\tDP\t234\n")
 
            file1 = os.path.join(input_dir.path, "A.snvs.vcf")
            file2 = os.path.join(input_dir.path, "A.indels.vcf")
            all_variants = merge_data([file1, file2])
 
            self.assertEqual('1\t2352\tA\tG\tfoo\tDP\t234\n', all_variants[0])
            self.assertEqual('1\t235234\tA\tG\tfoo\tDP\t234\n', all_variants[1])
            self.assertEqual('2\t2352\tA\tG\tfoo\tDP\t234\n', all_variants[2])
            self.assertEqual('1\t2700\tA\tG\tfoo\tDP\t345\n', all_variants[3])
            self.assertEqual('10\t2352\tA\tG\tfoo\tDP\t234\n', all_variants[4])
            self.assertEqual('1\t2\tA\tG\tfoo\tDP\t234\n', all_variants[5])
        input_dir.cleanup()
         
    def test_merge_mergeDataSamePos(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.snvs.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tG\tfoo\tDP\t234\n")
            input_dir.write("A.indels.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tGT\tfoo\tDP\t234\n")
 
            file1 = os.path.join(input_dir.path, "A.snvs.vcf")
            file2 = os.path.join(input_dir.path, "A.indels.vcf")
            all_variants = merge_data([file1, file2])
 
            self.assertEqual('1\t2352\tA\tG\tfoo\tDP\t234\n', all_variants[0])
            self.assertEqual('1\t2352\tA\tGT\tfoo\tDP\t234\n', all_variants[1])
        input_dir.cleanup()

class ValidateSplitLine(unittest.TestCase):
    def test_validateSplitLine_valid(self):
        split_line = ["chr1", "234", ".", "A", "G", ".", "PASS", "SS=5", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(0, invalid)
        self.assertEqual(0, warn)
        self.assertEqual(split_line, line)
         
        split_line = ["chr1", "234", ".", "A", "G", ".", "PASS", "SS=2", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(0, invalid)
        self.assertEqual(0, warn)
        self.assertEqual(split_line, line)
        
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

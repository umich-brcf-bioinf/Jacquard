from StringIO import StringIO
import sys
import os
import unittest
from testfixtures import TempDirectory
from jacquard.normalize import identify_merge_candidates, get_headers, merge_data, validate_split_line, _build_vcf_readers
from jacquard.variant_callers import varscan, strelka
    
class NormalizeTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr
        
    def test_indentifyMergeCandidates_Strelka(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/normalize_strelka_test/input/"
        in_files = [input_dir + "tiny_strelka.indels.vcf", input_dir + "tiny_strelka.snvs.vcf"]
        output_dir = script_dir + "/normalize_strelka_test/output/"

        merge_candidates, hc_candidates = identify_merge_candidates(in_files, output_dir, strelka.Strelka())
        
        self.assertEqual([output_dir + "tiny_strelka.merged.vcf"], merge_candidates.keys())
        self.assertEqual([[input_dir + "tiny_strelka.indels.vcf", input_dir + "tiny_strelka.snvs.vcf"]], merge_candidates.values())
        
    def test_indentifyMergeCandidates_missingFiles_VarScan(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/normalize_varscan_test/input/"
        in_files = [input_dir + "foo_indel.vcf", input_dir + "tiny_indel.vcf", input_dir + "tiny_snp.vcf", input_dir + "tiny_indel.Germline.hc", input_dir + "tiny_indel.LOH.hc", input_dir + "tiny_indel.Somatic.hc", input_dir + "tiny_snp.Germline.hc", input_dir + "tiny_snp.LOH.hc", input_dir + "tiny_snp.Somatic.hc"]
        output_dir = script_dir + "/normalize_varscan_test/output"
        
        with self.assertRaises(SystemExit) as cm:

            merge_candidates = identify_merge_candidates(in_files, output_dir, varscan.Varscan())
        self.assertEqual(cm.exception.code, 1)
        
    def test_indentifyMergeCandidates_VarScan(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/normalize_varscan_test/input/"
        in_files = [input_dir + "tiny_indel.vcf", input_dir + "tiny_snp.vcf", input_dir + "tiny_indel.Germline.hc", input_dir + "tiny_indel.LOH.hc", input_dir + "tiny_indel.Somatic.hc", input_dir + "tiny_snp.Germline.hc", input_dir + "tiny_snp.LOH.hc", input_dir + "tiny_snp.Somatic.hc"]
        output_dir = script_dir + "/normalize_varscan_test/output/"

        merge_candidates, hc_candidates = identify_merge_candidates(in_files, output_dir, varscan.Varscan())
        
        self.assertEqual([output_dir + "tiny_merged.vcf"], merge_candidates.keys())
        self.assertEqual([[input_dir + "tiny_indel.vcf", input_dir + "tiny_snp.vcf"]], merge_candidates.values())
    
    def test_build_vcf_readers(self):
        vcf_content ='''##source=strelka
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
'''
        vcf_content = vcf_content.replace('|',"\t")
        
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf",vcf_content)
            input_dir.write("B.vcf",vcf_content)
            
            vcf_readers = _build_vcf_readers(input_dir.path)
            
            self.assertEqual("A.vcf", vcf_readers[0].file_name)
            self.assertEqual(["##source=strelka"], vcf_readers[0].metaheaders)
            self.assertEqual("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR",
                              vcf_readers[0].column_header)
            self.assertEqual("B.vcf", vcf_readers[1].file_name)
            self.assertEqual(["##source=strelka"], vcf_readers[1].metaheaders)
            self.assertEqual("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR",
                              vcf_readers[1].column_header)
            self.assertEqual(2, len(vcf_readers))
        
    def test_merge_and_sort(self):
        pass
        
    def test_merge_getHeaders_Strelka(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf","##source=strelka\n##foobarbaz\n#CHROM\tNORMAL\tTUMOR\n123\n456\n")
 
            file = os.path.join(input_dir.path, "A.vcf")
            meta_headers, header = get_headers(file)
            self.assertEqual(["##source=strelka\n", "##foobarbaz\n"], meta_headers)
            self.assertEqual("#CHROM\tNORMAL\tTUMOR\n", header)
             
        input_dir.cleanup()
#     
    def test_merge_mergeData_Strelka(self):
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
         
    def test_merge_mergeDataSamePos_Strelka(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.snvs.vcf","##source=strelka\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tG\tfoo\tDP\t234\n")
            input_dir.write("A.indels.vcf","##source=strelka\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tGT\tfoo\tDP\t234\n")
 
            file1 = os.path.join(input_dir.path, "A.snvs.vcf")
            file2 = os.path.join(input_dir.path, "A.indels.vcf")
            all_variants = merge_data([file1, file2])
 
            self.assertEqual('1\t2352\tA\tG\tfoo\tDP\t234\n', all_variants[0])
            self.assertEqual('1\t2352\tA\tGT\tfoo\tDP\t234\n', all_variants[1])
        input_dir.cleanup()
        
    def test_merge_getHeaders_VarScan(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf","##source=VarScan2\n##foobarbaz\n#CHROM\tNORMAL\tTUMOR\n123\n456\n")

            file = os.path.join(input_dir.path, "A.vcf")
            meta_headers, header = get_headers(file)
            self.assertEqual(["##source=VarScan2\n", "##foobarbaz\n"], meta_headers)
            self.assertEqual("#CHROM\tNORMAL\tTUMOR\n", header)
            
        input_dir.cleanup()
    
    def test_merge_mergeData_VarScan(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.snp.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tG\tfoo\tDP\t234\n1\t235234\tA\tG\tfoo\tDP\t234\n2\t2352\tA\tG\tfoo\tDP\t234\n")
            input_dir.write("A.indel.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2700\tA\tG\tfoo\tDP\t345\n10\t2352\tA\tG\tfoo\tDP\t234\n1\t2\tA\tG\tfoo\tDP\t234\n")

            file1 = os.path.join(input_dir.path, "A.snp.vcf")
            file2 = os.path.join(input_dir.path, "A.indel.vcf")
            all_variants = merge_data([file1, file2])

            self.assertEqual('1\t2352\tA\tG\tfoo\tDP\t234\n', all_variants[0])
            self.assertEqual('1\t235234\tA\tG\tfoo\tDP\t234\n', all_variants[1])
            self.assertEqual('2\t2352\tA\tG\tfoo\tDP\t234\n', all_variants[2])
            self.assertEqual('1\t2700\tA\tG\tfoo\tDP\t345\n', all_variants[3])
            self.assertEqual('10\t2352\tA\tG\tfoo\tDP\t234\n', all_variants[4])
            self.assertEqual('1\t2\tA\tG\tfoo\tDP\t234\n', all_variants[5])
        input_dir.cleanup()
        
    def test_merge_mergeDataSamePos_VarScan(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.snp.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tG\tfoo\tDP\t234\n")
            input_dir.write("A.indel.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2352\tA\tGT\tfoo\tDP\t234\n")

            file1 = os.path.join(input_dir.path, "A.snp.vcf")
            file2 = os.path.join(input_dir.path, "A.indel.vcf")
            all_variants = merge_data([file1, file2])

            self.assertEqual('1\t2352\tA\tG\tfoo\tDP\t234\n', all_variants[0])
            self.assertEqual('1\t2352\tA\tGT\tfoo\tDP\t234\n', all_variants[1])
        input_dir.cleanup()

    def test_validateSplitLine_valid_Strelka(self):
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
        
    def test_validateSplitLine_valid_VarScan(self):
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
        
    def test_validateSplitLine_invalidAltOkaySS_VarScan(self):
        split_line = ["chr1", "234", ".", "A", "G/C", ".", "PASS", "SS=5", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(0, invalid)
        self.assertEqual(1, warn)
        self.assertEqual([], line)
        
        split_line = ["chr1", "234", ".", "A", "-G", ".", "PASS", "SS=5", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(0, invalid)
        self.assertEqual(1, warn)
        self.assertEqual([], line)
        
        split_line = ["chr1", "234", ".", "A", "+G", ".", "PASS", "SS=5", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(0, invalid)
        self.assertEqual(1, warn)
        self.assertEqual([], line)
        
    def test_validateSplitLine_invalidRefOkaySS_VarScan(self):
        split_line = ["chr1", "234", ".", "A/C", "G", ".", "PASS", "SS=5", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(0, invalid)
        self.assertEqual(1, warn)
        self.assertEqual([], line)
        
        split_line = ["chr1", "234", ".", "-A", "G", ".", "PASS", "SS=5", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(0, invalid)
        self.assertEqual(1, warn)
        self.assertEqual([], line)
        
        split_line= ["chr1", "234", ".", "+A", "G", ".", "PASS", "SS=5", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(0, invalid)
        self.assertEqual(1, warn)
        self.assertEqual([], line)
        
    def test_validateSplitLine_invalidAltBadSS_VarScan(self):
        split_line = ["chr1", "234", ".", "A", "G/C", ".", "PASS", "SS=2", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(1, invalid)
        self.assertEqual(0, warn)
        self.assertEqual([], line)
        
        split_line = ["chr1", "234", ".", "A", "+G", ".", "PASS", "SS=2", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(1, invalid)
        self.assertEqual(0, warn)
        self.assertEqual([], line)
        
        split_line = ["chr1", "234", ".", "A", "-G", ".", "PASS", "SS=2", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(1, invalid)
        self.assertEqual(0, warn)
        self.assertEqual([], line)
        
    def test_validateSplitLine_invalidRefBadSS_VarScan(self):
        split_line = ["chr1", "234", ".", "A/C", "G", ".", "PASS", "SS=2", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(1, invalid)
        self.assertEqual(0, warn)
        self.assertEqual([], line)
        
        split_line = ["chr1", "234", ".", "-A", "G", ".", "PASS", "SS=2", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(1, invalid)
        self.assertEqual(0, warn)
        self.assertEqual([], line)
        
        split_line = ["chr1", "234", ".", "+A", "G", ".", "PASS", "SS=2", "DP", "345"]
        invalid, warn, line = validate_split_line(split_line, 0, 0)
        self.assertEqual(1, invalid)
        self.assertEqual(0, warn)
        self.assertEqual([], line)
        
# class MockWriter():
#     def __init__(self):
#         self._content = []
#         self.wasClosed = False
# 
#     def write(self, content):
#         self._content.extend(content.splitlines())
#         
#     def lines(self):
#         return self._content
# 
#     def close(self):
#         self.wasClosed = True

#!/usr/bin/python2.7
from collections import defaultdict
import os
import unittest
import testfixtures
from testfixtures import TempDirectory
from bin.normalize_varscan import identify_merge_candidates, validate_directories, get_headers, merge_data, sort_data, change_pos_to_int, write_output, validate_split_line, identify_hc_variants, mark_hc_variants, validate_file_set

class IdentifyMergeCandidatesTestCase(unittest.TestCase):
    def test_indentifyMergeCandidates_missingFiles(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/normalize_varscan_test/input/"
        in_files = [input_dir + "foo_indel.vcf", input_dir + "tiny_indel.vcf", input_dir + "tiny_snp.vcf", input_dir + "tiny_indel.Germline.hc", input_dir + "tiny_indel.LOH.hc", input_dir + "tiny_indel.Somatic.hc", input_dir + "tiny_snp.Germline.hc", input_dir + "tiny_snp.LOH.hc", input_dir + "tiny_snp.Somatic.hc"]
        
        with self.assertRaises(SystemExit) as cm:
            merge_candidates = identify_merge_candidates(in_files)
        self.assertEqual(cm.exception.code, 1)
        
    def test_indentifyMergeCandidates_HC(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/normalize_varscan_test/input/"
        in_files = [input_dir + "tiny_indel.vcf", input_dir + "tiny_snp.vcf", input_dir + "tiny_indel.Germline.hc", input_dir + "tiny_indel.LOH.hc", input_dir + "tiny_indel.Somatic.hc", input_dir + "tiny_snp.Germline.hc", input_dir + "tiny_snp.LOH.hc", input_dir + "tiny_snp.Somatic.hc"]
        
        merge_candidates = identify_merge_candidates(in_files)
        
        self.assertEqual(["tiny_merged.vcf"], merge_candidates.keys())
        self.assertEqual([[input_dir + "tiny_indel.vcf", input_dir + "tiny_snp.vcf"]], merge_candidates.values())
        
    def test_identifyHcVariants(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/normalize_varscan_test/input/"
        hc_candidates = {"tiny_merged.vcf" : [input_dir + "tiny_indel.Germline.hc", input_dir + "tiny_indel.LOH.hc", input_dir + "tiny_indel.Somatic.hc", input_dir + "tiny_snp.Germline.hc", input_dir + "tiny_snp.LOH.hc", input_dir + "tiny_snp.Somatic.hc"]}
        
        hc_variants = identify_hc_variants(hc_candidates)
        expected_hc_variants = {'chr1^161332554^A^G': 1, 'chr9^33794812^G^T': 1, 'chr2^27015610^C^-CA': 1, 'chr3^156746013^T^C': 1, 'chr1^153773150^C^+T': 1, 'chr11^48387683^G^-T': 1, 
                                'chr11^1651198^G^-AGGCTGTGGGGGCTGTGGCTCCGGCTGTGC': 1, 'chr11^1651585^C^-CTGCTGCCAGTCCAGCTGCTGTAAGCCTTA': 1, 'chr2^32843743^A^+T': 1, 'chr9^33796627^T^C': 1, 
                                'chr2^33764141^G^+T': 1, 'chr1^16862212^C^T': 1, 'chr9^118950132^G^A': 1, 'chr1^2422614^G^A': 1, 'chr6^34004006^G^A': 1, 'chr5^139909358^A^C': 1, 
                                'chr1^180944532^T^+AA': 1, 'chr1^226044489^G^+C': 1, 'chr1^246939610^A^+TATT': 1, 'chr1^1575836^C^G': 1, 'chr2^27656822^G^+A': 1, 'chr2^31351842^A^+AC': 1, 
                                'chr1^158913532^C^T': 1, 'chr1^86487992^A^+ACAG': 1, 'chr9^33796630^T^C': 1, 'chr1^12856111^C^T': 1, 'chr6^36984880^C^-T': 1, 'chr2^25463303^G^-A': 1, 
                                'chr1^29446332^A^-T': 1, 'chr2^1133229^C^G': 1, 'chr1^7315580^C^T': 1, 'chr1^161332557^G^A': 1, 'chr1^171085314^A^+T': 1, 'chr1^204586805^G^A': 1, 
                                'chr2^24110911^G^+A': 1, 'chr1^14397^CTGT^C': 1, 'chr1^16890760^G^A': 1, 'chr2^242800^T^C': 1, 'chr1^172362744^A^+AAG': 1, 'chr1^12942306^C^T': 1, 'chr1^1886859^C^T': 1, 
                                'chr3^71026162^G^+T': 1, 'chr1^114137013^C^+AT': 1, 'chr2^19551276^T^-A': 1, 'chr1^762601^T^C': 1, 'chr1^149577615^A^G': 1, 'chr2^1442551^A^-T': 1, 'chr1^12898445^C^A': 1, 
                                'chr2^276942^A^G': 1, 'chr1^108160260^G^-A': 1, 'chr1^67292473^A^+T': 1, 'chr1^27661755^T^+GTA': 1, 'chr7^142498699^A^G': 1, 'chr2^99463179^G^A': 1, 'chr1^146013577^A^C': 1, 
                                'chr1^3638593^G^A': 1, 'chr9^18753463^G^A': 1, 'chr2^43778857^G^-A': 1, 'chr2^27015141^G^+C': 1, 'chr11^2428307^C^-CTCGGCCTCACCCAGGTGCTCCCGCTTGTG': 1, 'chr2^269352^G^A': 1, 
                                'chr1^38049373^C^T': 1, 'chr2^39025043^G^-A': 1, 'chr2^234130^T^G': 1, 'chr3^110837677^C^T': 1, 'chr1^15541607^T^C': 1, 'chr2^12877501^G^+A': 1, 'chr1^17020146^A^G': 1, 
                                'chr2^672745^T^C': 1, 'chr2^25061270^G^-GGGTGGGGT': 1, 'chr1^102462446^G^-TCAGT': 1, 'chr1^16914160^G^A': 1, 'chr2^231115^A^G': 1, 'chr1^45116470^C^T': 1, 'chr2^26678117^T^+A': 1, 
                                'chr17^79084047^C^-GGACCGCG': 1, 'chr9^33794809^G^A': 1, 'chr1^1647871^T^C': 1, 'chr9^33798170^C^A': 1, 'chr2^27170470^C^+G': 1, 'chr2^11682754^C^+G': 1, 'chr2^1093965^G^A': 1, 
                                'chr1^3607520^G^A': 1, 'chr1^172378884^G^+T': 1, 'chr1^86146672^G^A': 1, 'chr17^7190430^A^+G': 1, 'chr2^29158317^C^+T': 1, 'chr1^1647814^T^C': 1, 'chr7^45143090^G^+GACAGCC': 1, 
                                'chr1^16748087^G^A': 1, 'chr1^248801778^A^T': 1, 'chr2^243504^C^T': 1, 'chr1^51061718^T^-A': 1, 'chr7^91627046^C^G': 1, 'chr4^166964590^T^A': 1, 'chr4^88536899^A^G': 1, 'chr2^279705^C^T': 1}
        self.assertEqual(expected_hc_variants, hc_variants)
        
    def test_markHcCandidates(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/normalize_varscan_test/input/"
        hc_variants = {'chr1^161332554^A^G': 1, 'chr9^33794812^G^T': 1, 'chr2^27015610^C^-CA': 1, 'chr3^156746013^T^C': 1, 'chr1^153773150^C^+T': 1, 'chr11^48387683^G^-T': 1, 
                        'chr11^1651198^G^-AGGCTGTGGGGGCTGTGGCTCCGGCTGTGC': 1, 'chr11^1651585^C^-CTGCTGCCAGTCCAGCTGCTGTAAGCCTTA': 1, 'chr2^32843743^A^+T': 1, 'chr9^33796627^T^C': 1, 
                        'chr2^33764141^G^+T': 1, 'chr1^16862212^C^T': 1, 'chr9^118950132^G^A': 1, 'chr1^2422614^G^A': 1, 'chr6^34004006^G^A': 1, 'chr5^139909358^A^C': 1, 
                        'chr1^180944532^T^+AA': 1, 'chr1^226044489^G^+C': 1, 'chr1^246939610^A^+TATT': 1, 'chr1^1575836^C^G': 1, 'chr2^27656822^G^+A': 1, 'chr2^31351842^A^+AC': 1, 
                        'chr1^158913532^C^T': 1, 'chr1^86487992^A^+ACAG': 1, 'chr9^33796630^T^C': 1, 'chr1^12856111^C^T': 1, 'chr6^36984880^C^-T': 1, 'chr2^25463303^G^-A': 1, 
                        'chr1^29446332^A^-T': 1, 'chr2^1133229^C^G': 1, 'chr1^7315580^C^T': 1, 'chr1^161332557^G^A': 1, 'chr1^171085314^A^+T': 1, 'chr1^204586805^G^A': 1, 
                        'chr2^24110911^G^+A': 1, 'chr1^14397^CTGT^C': 1, 'chr1^16890760^G^A': 1, 'chr2^242800^T^C': 1, 'chr1^172362744^A^+AAG': 1, 'chr1^12942306^C^T': 1, 'chr1^1886859^C^T': 1, 
                        'chr3^71026162^G^+T': 1, 'chr1^114137013^C^+AT': 1, 'chr2^19551276^T^-A': 1, 'chr1^762601^T^C': 1, 'chr1^149577615^A^G': 1, 'chr2^1442551^A^-T': 1, 'chr1^12898445^C^A': 1, 
                        'chr2^276942^A^G': 1, 'chr1^108160260^G^-A': 1, 'chr1^67292473^A^+T': 1, 'chr1^27661755^T^+GTA': 1, 'chr7^142498699^A^G': 1, 'chr2^99463179^G^A': 1, 'chr1^146013577^A^C': 1, 
                        'chr1^3638593^G^A': 1, 'chr9^18753463^G^A': 1, 'chr2^43778857^G^-A': 1, 'chr2^27015141^G^+C': 1, 'chr11^2428307^C^-CTCGGCCTCACCCAGGTGCTCCCGCTTGTG': 1, 'chr2^269352^G^A': 1, 
                        'chr1^38049373^C^T': 1, 'chr2^39025043^G^-A': 1, 'chr2^234130^T^G': 1, 'chr3^110837677^C^T': 1, 'chr1^15541607^T^C': 1, 'chr2^12877501^G^+A': 1, 'chr1^17020146^A^G': 1, 
                        'chr2^672745^T^C': 1, 'chr2^25061270^G^-GGGTGGGGT': 1, 'chr1^102462446^G^-TCAGT': 1, 'chr1^16914160^G^A': 1, 'chr2^231115^A^G': 1, 'chr1^45116470^C^T': 1, 'chr2^26678117^T^+A': 1, 
                        'chr17^79084047^C^-GGACCGCG': 1, 'chr9^33794809^G^A': 1, 'chr1^1647871^T^C': 1, 'chr9^33798170^C^A': 1, 'chr2^27170470^C^+G': 1, 'chr2^11682754^C^+G': 1, 'chr2^1093965^G^A': 1, 
                        'chr1^3607520^G^A': 1, 'chr1^172378884^G^+T': 1, 'chr1^86146672^G^A': 1, 'chr17^7190430^A^+G': 1, 'chr2^29158317^C^+T': 1, 'chr1^1647814^T^C': 1, 'chr7^45143090^G^+GACAGCC': 1, 
                        'chr1^16748087^G^A': 1, 'chr1^248801778^A^T': 1, 'chr2^243504^C^T': 1, 'chr1^51061718^T^-A': 1, 'chr7^91627046^C^G': 1, 'chr4^166964590^T^A': 1, 'chr4^88536899^A^G': 1, 'chr2^279705^C^T': 1}
        merge_candidates = {"tiny_merged.vcf" : [input_dir + "tiny_snp.vcf", input_dir + "tiny_indel.vcf"]}
        
        actual_merge_candidates, marked_as_hc = mark_hc_variants(hc_variants, merge_candidates)
        
        self.assertEqual(merge_candidates, actual_merge_candidates)
        self.assertEqual(['chr1^14397^CTGT^C'], marked_as_hc)
        
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
        
    def test_merge_validateFileSet(self):
        all_keys = ["foo_merged.vcf", "foo_merged.Somatic.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
        sample_files = validate_file_set(all_keys)
        
        self.assertEqual({'foo': ['.vcf', '.Somatic.hc', '.Germline.hc', '.LOH.hc']}, sample_files)
        
    def test_merge_validateFileSet_invalid(self):
        all_keys = ["foo_merged.vcf", "foo_merged.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
        
        with self.assertRaises(SystemExit) as cm:
            sample_files = validate_file_set(all_keys)
        self.assertEqual(cm.exception.code, 1)
    
    def test_merge_validateFileSet_extraFiles(self):
        all_keys = ["foo_merged.bar.hc", "foo_merged.vcf", "foo_merged.Somatic.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
        sample_files = validate_file_set(all_keys)
        
        self.assertEqual({'foo': ['.bar.hc', '.vcf', '.Somatic.hc', '.Germline.hc', '.LOH.hc']}, sample_files)
        

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
        
    def test_validateSplitLine_invalidAltOkaySS(self):
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
        
    def test_validateSplitLine_invalidRefOkaySS(self):
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
        
    def test_validateSplitLine_invalidAltBadSS(self):
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
        
    def test_validateSplitLine_invalidRefBadSS(self):
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
            validate_directories(input_dir, first_out_dir + "/bar")
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

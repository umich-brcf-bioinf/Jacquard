# pylint: disable=C0103,C0301,R0904
from collections import OrderedDict,defaultdict
import os
import unittest

from jacquard.variant_callers import varscan
from jacquard.utils import __version__,jq_af_tag,jq_dp_tag,jq_somatic_tag,\
    JQException
from test.variant_callers.mutect_test import MockTag
from jacquard.vcf import VcfRecord 


class MockWriter():
    def __init__(self):
        self._content = []
        self.opened = False
        self.closed = False

    def open (self):
        self.opened = True
        
    def write(self, content):
        self._content.extend(content.splitlines())
        
    def lines(self):
        return self._content

    def close(self):
        self.closed = True
        
class MockFileReader(object):
    def __init__(self, input_filepath="/foo/mockFileReader.txt", content = []):
        self.input_filepath = input_filepath
        self.file_name = os.path.basename(input_filepath)
        self._content = content
        self.open_was_called = False
        self.close_was_called = False
    
    def open(self):
        self.open_was_called = True
    
    def read_lines(self):
        for line in self._content:
            yield line
         
    def close(self):
        self.close_was_called = True
        
        
class AlleleFreqTagTestCase(unittest.TestCase):

    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}VS,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ)",Source="Jacquard",Version={1}>'.format(jq_af_tag, __version__), varscan.AlleleFreqTag().metaheader)
                
    def test_format_missingAFTag(self):
        tag = varscan.AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        originalVcfRecord = VcfRecord(line)
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())
        
    def test_format_presentAFTag(self):
        tag = varscan.AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3|56.7%:SA.2:SA.3|83.4%:SB.2:SB.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3:JQ_AF_VS|56.7%:SA.2:SA.3:0.57|83.4%:SB.2:SB.3:0.83\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
    def test_format_multAlt(self):
        tag = varscan.AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3|56.7%,41%:SA.2:SA.3|83.4%,23%:SB.2:SB.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3:JQ_AF_VS|56.7%,41%:SA.2:SA.3:0.57,0.41|83.4%,23%:SB.2:SB.3:0.83,0.23\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
class DepthTagTestCase(unittest.TestCase):
    
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}VS,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP)",Source="Jacquard",Version={1}>'.format(jq_dp_tag, __version__), varscan.DepthTag().metaheader)
        
    def test_format_missingDPTag(self):
        tag = varscan.DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        originalVcfRecord = VcfRecord(line)
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())
        
    def test_format_presentDPTag(self):
        tag = varscan.DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3|2:SA.2:SA.3|4:SB.2:SB.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3:JQ_DP_VS|2:SA.2:SA.3:2|4:SB.2:SB.3:4\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
class SomaticTagTestCase(unittest.TestCase):
 
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1=somatic (based on SOMATIC info tag and if sample is TUMOR)",Source="Jacquard",Version={1}>'.format(jq_somatic_tag,__version__), varscan.SomaticTag().metaheader)
        
    def test_format_missingSSTag(self):
        tag = varscan.SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:"+jq_somatic_tag+"VS|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
    def test_format_presentSSTag_withHC(self):
        tag = varscan.SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|JQ_HC_VS;SS=2|F1:F2:F3|2:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|JQ_HC_VS;SS=2|F1:F2:F3:"+jq_somatic_tag+"VS|2:SA.2:SA.3:0|SB.1:SB.2:SB.3:1\n").replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
    def test_format_presentSSTag_withoutHC(self):
        tag = varscan.SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|SS=2|F1:F2:F3|2:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|SS=2|F1:F2:F3:"+jq_somatic_tag+"VS|2:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
    
    def test_format_SSTag_notEqual2(self):
        tag = varscan.SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|SS=5|F1:F2:F3|2:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|SS=5|F1:F2:F3:"+jq_somatic_tag+"VS|2:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class MockTag(object):
    def __init__(self, field_name, field_value, metaheader=None):
        self.field_name = field_name
        self.field_value = field_value
        self.metaheader = metaheader
    
    def format(self, vcfRecord):
        vcfRecord.insert_format_field(self.field_name, {0:self.field_value, 1:self.field_value})
                
class VarscanTestCase(unittest.TestCase):
     
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.caller = varscan.Varscan()
         
    def test_normalize(self):
        writer = MockWriter()
        content1 = ["##foo", "##bar", "#baz"]
        content2 = ["##hi", "##bar", "#baz"]
        readers = []
        readers.append(MockFileReader("indel.vcf", content1))
        readers.append(MockFileReader("snp.vcf", content2))
        self.append_hc_files(readers)
        self.caller.normalize(writer, readers)
         
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)
        self.assertEquals(["##bar", "##foo", "##hi", "#baz"], writer.lines())
        
    def test_normalize_wrongNumberOfFiles(self):
        self.assertRaisesRegexp(JQException,
            r"ERROR: VarScan directories should have exactly two input VCF files per patient, but found \[1\].",
            self.caller.normalize, MockWriter(), [MockFileReader(input_filepath="foo.vcf", content=["##foo", "#bar"])])
            
    def test_normalize_mismatchedColumnHeaders(self):
        writer = MockWriter()
        content1 = ["##foo", "##bar", "#baz"]
        content2 = ["##foo", "##bar", "#bluh"]
        readers = []
        readers.append(MockFileReader("indel.vcf", content1))
        readers.append(MockFileReader("snp.vcf", content2))
        self.append_hc_files(readers)
         
        self.assertRaisesRegexp(JQException, r"ERROR: The column headers for VCF files \[indel.vcf,snp.vcf\] do not match.",
                                 self.caller.normalize, writer, readers)
          
    def test_normalize_raisesExceptionMissingIndelSnvs(self):
        error_msg = r"ERROR: Each patient in a VarScan directory should have a snp file and an indel file."
        self.assert_two_files_throw_exception("foo.vcf", "bar.vcf", error_msg)
        self.assert_two_files_throw_exception("snp.vcf", "bar.vcf", error_msg)
        self.assert_two_files_throw_exception("foo.snp.vcf", "bar.vcf", error_msg)
        self.assert_two_files_throw_exception("foo.indel.vcf", "bar.vcf", error_msg)
        self.assert_two_files_throw_exception("foo.indel.vcf", "bar.indel.vcf", error_msg)
        self.assert_two_files_throw_exception("foo.snp.vcf", "bar.snp.vcf", error_msg)
        self.assert_two_files_throw_exception("snp/foo.vcf", "indel/bar.vcf", error_msg)
        self.assert_two_files_throw_exception("indel.snp.vcf", "bar.vcf", error_msg)
        self.assert_two_files_throw_exception("A.indel.snp.vcf", "B.indel.snp.vcf", error_msg)
         
    def test_normalize_writesSequentialRecords(self):
        writer = MockWriter()
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        record2 = "chr2\t.\t.\t.\t.\t.\t.\t.\t."
        record3 = "chr3\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "#bar", record2, record3]
        content2 = ["##foo", "#bar", record1, record3]
        readers = []
        readers.append(MockFileReader("indel.vcf", content1))
        readers.append(MockFileReader("snp.vcf", content2))
        self.append_hc_files(readers)
        self.caller.normalize(writer,readers)
           
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)
        self.assertEquals(["##foo", "#bar", record1, record2, record3, record3], writer.lines())
#         
    def assert_two_files_throw_exception(self, file1, file2, exception):
        readers = []
        content = ["##foo","#bar"]
        readers.append(MockFileReader(input_filepath=file1, content=content))
        readers.append(MockFileReader(input_filepath=file2, content=content))
        self.append_hc_files(readers)
    
#         self.assertRaises(JQException, self.caller.normalize, MockWriter(), readers)
        with self.assertRaisesRegexp(JQException,exception):
            self.caller.normalize(MockWriter(), readers)

    def append_hc_files(self,readers):
        readers.append(MockFileReader("snp.somatic.hc", []))
        readers.append(MockFileReader("indel.somatic.hc", []))
        
#     def setUp(self):
#         self.caller = varscan.Varscan()
#         
#     def test_validateInputFile_valid(self):
#         metaheaders = ["##source=VarScan2"]
#         valid = self.caller.validate_input_file(metaheaders,"#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR".replace('|',"\t"))
#          
#         self.assertEquals(1, valid)
#          
#     def test_validateInputFile_invalid_columnHeader(self):
#         metaheaders = ["##source=VarScan2"]
#         self.assertRaises(JQException, self.caller.validate_input_file, metaheaders,"#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|FOO|BAR".replace('|',"\t"))        
#     
#     def test_validateInputFile_invalid_metaheaders(self):
#         metaheaders = ["##foo"]
#         valid = self.caller.validate_input_file(metaheaders,"#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR".replace('|',"\t"))
#         self.assertEquals(0, valid)
#         
#     def test_finalSteps(self):
#         script_dir = os.path.dirname(os.path.abspath(__file__))
#         hc_candidates  = {os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_merged.Somatic.hc'): [
#             os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_indel.Somatic.hc'), os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_snp.Somatic.hc')], os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_merged.LOH.hc'): [os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_indel.LOH.hc'), os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_snp.LOH.hc')], os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_merged.Germline.hc'): [os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_indel.Germline.hc'), os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_snp.Germline.hc')]}
#         merge_candidates = {os.path.join(script_dir, '../reference_files/normalize_varscan_test/tiny_merged.vcf'): [os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_indel.vcf'), os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_snp.vcf')]}
#         output_dir = os.path.join(script_dir, "../reference_files/normalize_varscan_test/output")
#         marked_as_hc = self.caller.final_steps(hc_candidates, merge_candidates, output_dir)
#           
#         self.assertEquals(['chr1^14397^CTGT^C'], marked_as_hc)
#           
#     def test_handleHCFiles(self):
#         script_dir = os.path.dirname(os.path.abspath(__file__))
#         in_file = os.path.join(script_dir, "normalize_varscan_test", "input", "tiny_indel.Germline.hc")
#         out_dir = os.path.join(script_dir, "normalize_varscan_test")
#         hc_candidates = defaultdict(list)
#         hc_candidates = self.caller.handle_hc_files(in_file, out_dir, hc_candidates)
#           
#         expected_hc_candidates = {os.path.join(script_dir, "normalize_varscan_test", "tiny_merged.Germline.hc"): [os.path.join(script_dir, "normalize_varscan_test", "input", "tiny_indel.Germline.hc")]}
#         self.assertEquals(expected_hc_candidates, hc_candidates)
#               
#     def test_validateFileSet(self):
#         all_keys = ["foo_merged.vcf", "foo_merged.Somatic.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
#         sample_files = self.caller.validate_file_set(all_keys)
#           
#         self.assertEqual({'foo': ['.vcf', '.Somatic.hc', '.Germline.hc', '.LOH.hc']}, sample_files)
#           
#     def test_validateFileSet_invalid(self):
#         all_keys = ["foo_merged.vcf", "foo_merged.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
#           
#         with self.assertRaises(SystemExit) as cm:
#             self.caller.validate_file_set(all_keys)
#         self.assertEqual(cm.exception.code, 1)
#       
#     def test_validateFileSet_extraFiles(self):
#         all_keys = ["foo_merged.bar.hc", "foo_merged.vcf", "foo_merged.Somatic.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
#         sample_files = self.caller.validate_file_set(all_keys)
#           
#         self.assertEqual({'foo': ['.bar.hc', '.vcf', '.Somatic.hc', '.Germline.hc', '.LOH.hc']}, sample_files)
#          
#     def test_identifyHcVariants_VarScan(self):
#         script_dir = os.path.dirname(os.path.abspath(__file__))
#         input_dir = script_dir + "/../reference_files/normalize_varscan_test/input/"
#         hc_candidates = {"tiny_merged.vcf" : [input_dir + "tiny_indel.Germline.hc", input_dir + "tiny_indel.LOH.hc", input_dir + "tiny_indel.Somatic.hc", input_dir + "tiny_snp.Germline.hc", input_dir + "tiny_snp.LOH.hc", input_dir + "tiny_snp.Somatic.hc"]}
#            
#         hc_variants = self.caller.identify_hc_variants(hc_candidates)
#         expected_hc_variants = {'chr1^161332554^A^G': 1, 'chr9^33794812^G^T': 1, 'chr2^27015610^C^-CA': 1, 'chr3^156746013^T^C': 1, 'chr1^153773150^C^+T': 1, 'chr11^48387683^G^-T': 1, 
#                                 'chr11^1651198^G^-AGGCTGTGGGGGCTGTGGCTCCGGCTGTGC': 1, 'chr11^1651585^C^-CTGCTGCCAGTCCAGCTGCTGTAAGCCTTA': 1, 'chr2^32843743^A^+T': 1, 'chr9^33796627^T^C': 1, 
#                                 'chr2^33764141^G^+T': 1, 'chr1^16862212^C^T': 1, 'chr9^118950132^G^A': 1, 'chr1^2422614^G^A': 1, 'chr6^34004006^G^A': 1, 'chr5^139909358^A^C': 1, 
#                                 'chr1^180944532^T^+AA': 1, 'chr1^226044489^G^+C': 1, 'chr1^246939610^A^+TATT': 1, 'chr1^1575836^C^G': 1, 'chr2^27656822^G^+A': 1, 'chr2^31351842^A^+AC': 1, 
#                                 'chr1^158913532^C^T': 1, 'chr1^86487992^A^+ACAG': 1, 'chr9^33796630^T^C': 1, 'chr1^12856111^C^T': 1, 'chr6^36984880^C^-T': 1, 'chr2^25463303^G^-A': 1, 
#                                 'chr1^29446332^A^-T': 1, 'chr2^1133229^C^G': 1, 'chr1^7315580^C^T': 1, 'chr1^161332557^G^A': 1, 'chr1^171085314^A^+T': 1, 'chr1^204586805^G^A': 1, 
#                                 'chr2^24110911^G^+A': 1, 'chr1^14397^CTGT^C': 1, 'chr1^16890760^G^A': 1, 'chr2^242800^T^C': 1, 'chr1^172362744^A^+AAG': 1, 'chr1^12942306^C^T': 1, 'chr1^1886859^C^T': 1, 
#                                 'chr3^71026162^G^+T': 1, 'chr1^114137013^C^+AT': 1, 'chr2^19551276^T^-A': 1, 'chr1^762601^T^C': 1, 'chr1^149577615^A^G': 1, 'chr2^1442551^A^-T': 1, 'chr1^12898445^C^A': 1, 
#                                 'chr2^276942^A^G': 1, 'chr1^108160260^G^-A': 1, 'chr1^67292473^A^+T': 1, 'chr1^27661755^T^+GTA': 1, 'chr7^142498699^A^G': 1, 'chr2^99463179^G^A': 1, 'chr1^146013577^A^C': 1, 
#                                 'chr1^3638593^G^A': 1, 'chr9^18753463^G^A': 1, 'chr2^43778857^G^-A': 1, 'chr2^27015141^G^+C': 1, 'chr11^2428307^C^-CTCGGCCTCACCCAGGTGCTCCCGCTTGTG': 1, 'chr2^269352^G^A': 1, 
#                                 'chr1^38049373^C^T': 1, 'chr2^39025043^G^-A': 1, 'chr2^234130^T^G': 1, 'chr3^110837677^C^T': 1, 'chr1^15541607^T^C': 1, 'chr2^12877501^G^+A': 1, 'chr1^17020146^A^G': 1, 
#                                 'chr2^672745^T^C': 1, 'chr2^25061270^G^-GGGTGGGGT': 1, 'chr1^102462446^G^-TCAGT': 1, 'chr1^16914160^G^A': 1, 'chr2^231115^A^G': 1, 'chr1^45116470^C^T': 1, 'chr2^26678117^T^+A': 1, 
#                                 'chr17^79084047^C^-GGACCGCG': 1, 'chr9^33794809^G^A': 1, 'chr1^1647871^T^C': 1, 'chr9^33798170^C^A': 1, 'chr2^27170470^C^+G': 1, 'chr2^11682754^C^+G': 1, 'chr2^1093965^G^A': 1, 
#                                 'chr1^3607520^G^A': 1, 'chr1^172378884^G^+T': 1, 'chr1^86146672^G^A': 1, 'chr17^7190430^A^+G': 1, 'chr2^29158317^C^+T': 1, 'chr1^1647814^T^C': 1, 'chr7^45143090^G^+GACAGCC': 1, 
#                                 'chr1^16748087^G^A': 1, 'chr1^248801778^A^T': 1, 'chr2^243504^C^T': 1, 'chr1^51061718^T^-A': 1, 'chr7^91627046^C^G': 1, 'chr4^166964590^T^A': 1, 'chr4^88536899^A^G': 1, 'chr2^279705^C^T': 1}
#         self.assertEqual(expected_hc_variants, hc_variants)
#            
#     def test_markHcCandidates_VarScan(self):
#         script_dir = os.path.dirname(os.path.abspath(__file__))
#         input_dir = script_dir + "/../reference_files/normalize_varscan_test/merged_vcf/"
#         output_dir = script_dir + "/../reference_files/normalize_varscan_test/output/"
#         hc_variants = {'chr1^161332554^A^G': 1, 'chr9^33794812^G^T': 1, 'chr2^27015610^C^-CA': 1, 'chr3^156746013^T^C': 1, 'chr1^153773150^C^+T': 1, 'chr11^48387683^G^-T': 1, 
#                         'chr11^1651198^G^-AGGCTGTGGGGGCTGTGGCTCCGGCTGTGC': 1, 'chr11^1651585^C^-CTGCTGCCAGTCCAGCTGCTGTAAGCCTTA': 1, 'chr2^32843743^A^+T': 1, 'chr9^33796627^T^C': 1, 
#                         'chr2^33764141^G^+T': 1, 'chr1^16862212^C^T': 1, 'chr9^118950132^G^A': 1, 'chr1^2422614^G^A': 1, 'chr6^34004006^G^A': 1, 'chr5^139909358^A^C': 1, 
#                         'chr1^180944532^T^+AA': 1, 'chr1^226044489^G^+C': 1, 'chr1^246939610^A^+TATT': 1, 'chr1^1575836^C^G': 1, 'chr2^27656822^G^+A': 1, 'chr2^31351842^A^+AC': 1, 
#                         'chr1^158913532^C^T': 1, 'chr1^86487992^A^+ACAG': 1, 'chr9^33796630^T^C': 1, 'chr1^12856111^C^T': 1, 'chr6^36984880^C^-T': 1, 'chr2^25463303^G^-A': 1, 
#                         'chr1^29446332^A^-T': 1, 'chr2^1133229^C^G': 1, 'chr1^7315580^C^T': 1, 'chr1^161332557^G^A': 1, 'chr1^171085314^A^+T': 1, 'chr1^204586805^G^A': 1, 
#                         'chr2^24110911^G^+A': 1, 'chr1^14397^CTGT^C': 1, 'chr1^16890760^G^A': 1, 'chr2^242800^T^C': 1, 'chr1^172362744^A^+AAG': 1, 'chr1^12942306^C^T': 1, 'chr1^1886859^C^T': 1, 
#                         'chr3^71026162^G^+T': 1, 'chr1^114137013^C^+AT': 1, 'chr2^19551276^T^-A': 1, 'chr1^762601^T^C': 1, 'chr1^149577615^A^G': 1, 'chr2^1442551^A^-T': 1, 'chr1^12898445^C^A': 1, 
#                         'chr2^276942^A^G': 1, 'chr1^108160260^G^-A': 1, 'chr1^67292473^A^+T': 1, 'chr1^27661755^T^+GTA': 1, 'chr7^142498699^A^G': 1, 'chr2^99463179^G^A': 1, 'chr1^146013577^A^C': 1, 
#                         'chr1^3638593^G^A': 1, 'chr9^18753463^G^A': 1, 'chr2^43778857^G^-A': 1, 'chr2^27015141^G^+C': 1, 'chr11^2428307^C^-CTCGGCCTCACCCAGGTGCTCCCGCTTGTG': 1, 'chr2^269352^G^A': 1, 
#                         'chr1^38049373^C^T': 1, 'chr2^39025043^G^-A': 1, 'chr2^234130^T^G': 1, 'chr3^110837677^C^T': 1, 'chr1^15541607^T^C': 1, 'chr2^12877501^G^+A': 1, 'chr1^17020146^A^G': 1, 
#                         'chr2^672745^T^C': 1, 'chr2^25061270^G^-GGGTGGGGT': 1, 'chr1^102462446^G^-TCAGT': 1, 'chr1^16914160^G^A': 1, 'chr2^231115^A^G': 1, 'chr1^45116470^C^T': 1, 'chr2^26678117^T^+A': 1, 
#                         'chr17^79084047^C^-GGACCGCG': 1, 'chr9^33794809^G^A': 1, 'chr1^1647871^T^C': 1, 'chr9^33798170^C^A': 1, 'chr2^27170470^C^+G': 1, 'chr2^11682754^C^+G': 1, 'chr2^1093965^G^A': 1, 
#                         'chr1^3607520^G^A': 1, 'chr1^172378884^G^+T': 1, 'chr1^86146672^G^A': 1, 'chr17^7190430^A^+G': 1, 'chr2^29158317^C^+T': 1, 'chr1^1647814^T^C': 1, 'chr7^45143090^G^+GACAGCC': 1, 
#                         'chr1^16748087^G^A': 1, 'chr1^248801778^A^T': 1, 'chr2^243504^C^T': 1, 'chr1^51061718^T^-A': 1, 'chr7^91627046^C^G': 1, 'chr4^166964590^T^A': 1, 'chr4^88536899^A^G': 1, 'chr2^279705^C^T': 1}
#         merge_candidates = {input_dir + "tiny_merged.vcf" : [input_dir + "tiny_snp.vcf", input_dir + "tiny_indel.vcf"]}
#            
#         marked_as_hc = self.caller.mark_hc_variants(hc_variants, merge_candidates, output_dir)
#            
#         self.assertEqual(['chr1^14397^CTGT^C'], marked_as_hc)
#           

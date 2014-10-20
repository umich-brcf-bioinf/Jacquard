# pylint: disable=R0904,C0103
# pylint: disable=C0301

import unittest
import jacquard.variant_callers.strelka as strelka
from jacquard.utils import __version__, jq_af_tag, jq_somatic_tag, jq_dp_tag, JQException
import os

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
        self.assertEqual('##FORMAT=<ID={0}SK,Number=A,Type=Float,Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth. Uses TAR if available, otherwise uses uses DP2 if available, otherwise uses ACGT tier2 depth)",Source="Jacquard",Version={1}>'.format(jq_af_tag, __version__), strelka._AlleleFreqTag().metaheader)
     
    def test_format_missingAFTag(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        originalVcfRecord = VcfRecord(line)
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())
          
    def test_format_AUTag(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|A,C|QUAL|FILTER|INFO|AU:CU:GU:TU|1,2:3,4:5,6:7,8|9,10:11,12:13,14:15,16\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|A,C|QUAL|FILTER|INFO|AU:CU:GU:TU:JQ_AF_SK|1,2:3,4:5,6:7,8:0.1,0.2|9,10:11,12:13,14:15,16:0.19,0.23\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
    def test_format_AFTag_noAlt(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|.|QUAL|FILTER|INFO|AU:CU:GU:TU|1,2:3,4:5,6:7,8|9,10:11,12:13,14:15,16\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|.|QUAL|FILTER|INFO|AU:CU:GU:TU:JQ_AF_SK|1,2:3,4:5,6:7,8:.|9,10:11,12:13,14:15,16:.\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
    def test_format_TARTag(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:TAR|10:3,4|20:11,7\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:TAR:JQ_AF_SK|10:3,4:0.4|20:11,7:0.35\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
#         
class DepthTagTestCase(unittest.TestCase):
     
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}SK,Number=1,Type=Float,Description="Jacquard depth for Strelka (uses DP2 if available, otherwise uses ACGT tier2 depth)",Source="Jacquard",Version={1}>'.format(jq_dp_tag, __version__), strelka._DepthTag().metaheader)

    def test_format_missingTag(self):
        tag = strelka._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        originalVcfRecord = VcfRecord(line)
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())
         
    def test_format_DP2Tag(self):
        tag = strelka._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:F2:F3|2:SA.2:SA.3|4:SB.2:SB.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:F2:F3:JQ_DP_SK|2:SA.2:SA.3:2|4:SB.2:SB.3:4\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
    
    def test_format_AUTag(self):
        tag = strelka._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|AU:CU:GU:TU|1,2:3,4:5,6:7,8|9,10:11,12:13,14:15,16\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|AU:CU:GU:TU:JQ_DP_SK|1,2:3,4:5,6:7,8:20|9,10:11,12:13,14:15,16:52\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
         
class SomaticTagTestCase(unittest.TestCase):
     
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}SK,Number=1,Type=Integer,Description="Jacquard somatic status for Strelka: 0=non-somatic,1=somatic (based on PASS in FILTER column)",Source="Jacquard",Version={1}>'.format(jq_somatic_tag, __version__), strelka._SomaticTag().metaheader)
        
    def test_format_missingPASS(self):
        tag = strelka._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:"+jq_somatic_tag+"SK|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
         
    def test_format_presentPASS(self):
        tag = strelka._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|PASS|INFO|SS:F2:F3|2:SA.2:SA.3|5:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|PASS|INFO|SS:F2:F3:"+jq_somatic_tag+"SK|2:SA.2:SA.3:0|5:SB.2:SB.3:1\n").replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

# class MockTag(object):
#     def __init__(self, field_name, field_value, metaheader=None):
#         self.field_name = field_name
#         self.field_value = field_value
#         self.metaheader = metaheader
#     
#     def format(self, vcfRecord):
#         vcfRecord.insert_format_field(self.field_name, {0:self.field_value, 1:self.field_value})
#         
# class Mutect_TestCase(unittest.TestCase):
#     
#     def setUp(self):
#         self.caller = strelka.Mutect()
#         
#     def test_validateInputFile_isValid(self):
#         line = ["##MuTect"]
#         self.assertTrue(self.caller.validate_input_file(line))
#     
#     def test_validateInputFile_isNotValid(self):
#         line = ["Foo"]
#         self.assertFalse(self.caller.validate_input_file(line))
#     
#     def test_addTags(self):
#         input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
#         self.caller.tags=[MockTag("mockTag", 42)]
#         actual_line = self.caller.add_tags(VcfRecord(input_line))
# 
#         expected_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:mockTag|SA.1:SA.2:SA.3:42|SB.1:SB.2:SB.3:42\n".replace('|',"\t")
# 
#         self.assertEquals(expected_line, actual_line)
#         
#     def test_updateMetaheader(self):
#         input_metaheader = "#foo\n"
#         self.caller.tags=[MockTag("mockTag", 42, "##my_metaheader\n")]
#         actual_metaheader = self.caller.update_metaheader(input_metaheader)
# 
#         expected_line = "#foo\n##my_metaheader\n"
# 
#         self.assertEquals(expected_line, actual_metaheader)



class StrelkaTestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.caller = strelka.Strelka()
         
    def test_normalize(self):
        writer = MockWriter()
        content1 = ["##foo", "##bar", "#baz"]
        content2 = ["##hi", "##bar", "#baz"]
        reader1 = MockFileReader("indels.vcf", content1)
        reader2 = MockFileReader("snvs.vcf", content2)
        self.caller.normalize(writer,[reader1,reader2])
         
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)
        self.assertEquals(["##bar", "##foo", "##hi", "#baz"], writer.lines())
        
    def test_normalize_mismatchedColumnHeaders(self):
        writer = MockWriter()
        content1 = ["##foo", "##bar", "#baz"]
        content2 = ["##foo", "##bar", "#bluh"]
        reader1 = MockFileReader("indels.vcf", content1)
        reader2 = MockFileReader("snvs.vcf", content2)
        
        self.assertRaisesRegexp(JQException, r"The column headers for VCF files \[indels.vcf,snvs.vcf\] do not match.",
                                 self.caller.normalize, writer, [reader1,reader2])
        
#     def test_normalize_hasIndelSnvs(self):
#         writer = MockWriter()
# 
#         reader1 = MockFileReader("indels",["##metaheader\n","#column_header\n"])
#         reader2 = MockFileReader("snvs",["##metaheader\n","#column_header\n"])
#         
#         self.caller.normalize(writer,[reader1,reader2])
#         pass
    def test_normalize_wrongNumberOfFiles(self):
        self.assertRaisesRegexp(JQException,
                                r"Strelka directories should have exactly two input files per patient, but found \[1\].",
                                self.caller.normalize, MockWriter(), [MockFileReader(input_filepath="foo")])
            
    def test_normalize_raisesExceptionMissingIndelSnvs(self):
        self.assert_two_vcf_files_throw_exception("foo", "bar")
        self.assert_two_vcf_files_throw_exception("snvs", "bar")
        self.assert_two_vcf_files_throw_exception("foo.snvs", "bar")
        self.assert_two_vcf_files_throw_exception("foo.indels", "bar")
        self.assert_two_vcf_files_throw_exception("foo.indels", "bar.indels")
        self.assert_two_vcf_files_throw_exception("foo.snvs", "bar.snvs")
        self.assert_two_vcf_files_throw_exception("snvs/foo", "indels/bar")
        self.assert_two_vcf_files_throw_exception("indels.snvs", "bar")
        self.assert_two_vcf_files_throw_exception("A.indels.snvs", "B.indels.snvs")
        
    def test_normalize_writesSequentialRecords(self):
        writer = MockWriter()
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        record2 = "chr2\t.\t.\t.\t.\t.\t.\t.\t."
        record3 = "chr3\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "#bar", record2, record3]
        content2 = ["##foo", "#bar", record1, record3]
        reader1 = MockFileReader("indels.vcf", content1)
        reader2 = MockFileReader("snvs.vcf", content2)
        self.caller.normalize(writer,[reader1,reader2])
         
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)
        self.assertEquals(["##foo", "#bar", record1, record2, record3, record3], writer.lines())
        
    def assert_two_vcf_files_throw_exception(self, file1, file2):
        with self.assertRaisesRegexp(JQException,r"Each patient in a Strelka directory should have a snvs file and an indels file."):
            self.caller.normalize(MockWriter(), [MockFileReader(input_filepath=file1), MockFileReader(input_filepath=file2)])
            
#     def test_validateInputFile_valid(self):
#         caller = strelka.Strelka()
#         input_file = ["##source=strelka", "foo"]
#         name, valid = caller.validate_input_file(input_file)
# 
#         self.assertEquals("Strelka", name)
#         self.assertEquals(1, valid)
# 
#     def test_validateInputFile_invalid(self):
#         caller = strelka.Strelka()
#         input_file = ["##bar", "foo"]
#         name, valid = caller.validate_input_file(input_file)
# 
#         self.assertEquals("Strelka", name)
#         self.assertEquals(0, valid)
# 
#     def test_finalSteps(self):
#         caller = strelka.Strelka()
#         script_dir = os.path.dirname(os.path.abspath(__file__))
#         hc_candidates = {os.path.join(script_dir, 'normalize_varscan_test/input/tiny_merged.Somatic.hc'): [os.path.join(script_dir, 'normalize_varscan_test/input/tiny_indel.Somatic.hc'), os.path.join(script_dir, 'normalize_varscan_test/input/tiny_snp.Somatic.hc')], os.path.join(script_dir, 'normalize_varscan_test/input/tiny_merged.LOH.hc'): [os.path.join(script_dir, 'normalize_varscan_test/input/tiny_indel.LOH.hc'), os.path.join(script_dir, 'normalize_varscan_test/input/tiny_snp.LOH.hc')], os.path.join(script_dir, 'normalize_varscan_test/input/tiny_merged.Germline.hc'): [os.path.join(script_dir, 'normalize_varscan_test/input/tiny_indel.Germline.hc'), os.path.join(script_dir, 'normalize_varscan_test/input/tiny_snp.Germline.hc')]}
#         merge_candidates = {os.path.join(script_dir, 'normalize_varscan_test/output/tiny_merged.vcf'): [os.path.join(script_dir, 'normalize_varscan_test/input/tiny_indel.vcf'), os.path.join(script_dir, 'normalize_varscan_test/input/tiny_snp.vcf')]}
#         output_dir = os.path.join(script_dir, "normalize_varscan_test/output")
#         actual_merge_candidates = caller.final_steps(hc_candidates, merge_candidates, output_dir)
# 
#         self.assertEquals(merge_candidates, actual_merge_candidates)
# 
#     def test_handleHCFiles(self):
#         caller = strelka.Strelka()
#         script_dir = os.path.dirname(os.path.abspath(__file__))
#         in_file = os.path.join(script_dir, "normalize_varscan_test/input/tiny_indel.Germline.hc")
#         out_dir = os.path.join(script_dir, "normalize_varscan_test/output")
#         hc_candidates = defaultdict(list)
#         hc_candidates = caller.handle_hc_files(in_file, out_dir, hc_candidates)
# 
#         self.assertEquals(defaultdict(list), hc_candidates)
# 
#     def test_validateFileSet(self):
#         caller = strelka.Strelka()
#         all_keys = ["foo_merged.vcf", "foo_merged.Somatic.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
#         sample_files = caller.validate_file_set(all_keys)
# 
#         self.assertEqual(None, sample_files)

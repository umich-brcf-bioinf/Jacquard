# pylint: disable=R0904,C0103
# pylint: disable=C0301
from collections import OrderedDict, defaultdict
import os
import unittest

import jacquard.variant_callers.strelka as strelka
from jacquard.jacquard_utils import __version__, jq_af_tag
from jacquard.vcf_record import VcfRecord 

class AlleleFreqTagTestCase(unittest.TestCase):

    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}SK,Number=A,Type=Float,Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth. Uses TAR if available, otherwise uses uses DP2 if available, otherwise uses ACGT tier2 depth)",Source="Jacquard",Version={1}>\n'.format(jq_af_tag, __version__), strelka._AlleleFreqTag().metaheader)
    
#     def test_format_missingAFTag(self):
#         tag = strelka._AlleleFreqTag()
#         line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
#         originalVcfRecord = VcfRecord(line)
#         processedVcfRecord = VcfRecord(line)
#         tag.format(processedVcfRecord)
#         self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())
#          
#     def test_format_presentAFTag(self):
#         tag = strelka._AlleleFreqTag()
#         line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3|0.567:SA.2:SA.3|0.834:SB.2:SB.3\n".replace('|',"\t")
#         expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3:JQ_AF_MT|0.567:SA.2:SA.3:0.57|0.834:SB.2:SB.3:0.83\n".replace('|',"\t")
#         processedVcfRecord = VcfRecord(line)
#         tag.format(processedVcfRecord)
#         self.assertEquals(expected, processedVcfRecord.asText())
#         
#     def test_format_multAlt(self):
#         tag = strelka._AlleleFreqTag()
#         line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3|0.5,0.8:SA.2:SA.3|0.7,0.6:SB.2:SB.3\n".replace('|',"\t")
#         expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3:JQ_AF_MT|0.5,0.8:SA.2:SA.3:0.5,0.8|0.7,0.6:SB.2:SB.3:0.7,0.6\n".replace('|',"\t")
#         processedVcfRecord = VcfRecord(line)
#         tag.format(processedVcfRecord)
#         self.assertEquals(expected, processedVcfRecord.asText())
#         
class DepthTagTestCase(unittest.TestCase):
     
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
         
# class SomaticTagTestCase(unittest.TestCase):
#     
#     def test_format_missingSSTag(self):
#         tag = strelka._SomaticTag()
#         line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
#         expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:JQ_HC_SOM_MT|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n".replace('|',"\t")
#         processedVcfRecord = VcfRecord(line)
#         tag.format(processedVcfRecord)
#         self.assertEquals(expected, processedVcfRecord.asText())
#         
#     def test_format_presentSSTag(self):
#         tag = strelka._SomaticTag()
#         line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|SS:F2:F3|2:SA.2:SA.3|5:SB.2:SB.3\n".replace('|',"\t")
#         expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|SS:F2:F3:JQ_HC_SOM_MT|2:SA.2:SA.3:1|5:SB.2:SB.3:0\n".replace('|',"\t")
#         processedVcfRecord = VcfRecord(line)
#         tag.format(processedVcfRecord)
#         self.assertEquals(expected, processedVcfRecord.asText())
#         
#         
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


# class StrelkaTestCase(unittest.TestCase):
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
# 
# class Strelka_AlleleFreqTagTestCase(unittest.TestCase):
#     def test_metaheader(self):
#         self.assertEqual('##FORMAT=<ID=JQ_AF_SK,Number=A,Type=Float,Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth. Uses TAR if available, otherwise uses uses DP2 if available, otherwise uses ACGT tier2 depth)",Source="Jacquard",Version={0}>\n'.format(__version__), strelka.AlleleFreqTag().metaheader)
# 
#     def test_format_missingAFTag(self):
#         tag = strelka.AlleleFreqTag()
#         format_param_string = "A:B"
#         format_value_string = "1:2"
#         format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("C", "filter", "", format_dict, 0))
# 
#     def test_format_rounds(self):
#         tag = strelka.AlleleFreqTag()
#         format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:0.2,0.5:0.2,0.5:0.2,0.5:0.2,0.5:".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('AU', '0.2,0.5'), ('CU', '0.2,0.5'), ('GU', '0.2,0.5'), ('TU', '0.2,0.5'), ('JQ_AF_SK', '0.25')]), tag.format("C", "filter", "", format_dict, 0))
# 
#         format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:0.2,0.6:0.2,0.6:0.2,0.5:0.2,0.5:".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('AU', '0.2,0.6'), ('CU', '0.2,0.6'), ('GU', '0.2,0.5'), ('TU', '0.2,0.5'), ('JQ_AF_SK', '0.27')]), tag.format("C", "filter", "", format_dict, 0))
# 
#         format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:0.2,0.6:0.2,0.8:0.2,0.5:0.2,0.5:".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('AU', '0.2,0.6'), ('CU', '0.2,0.8'), ('GU', '0.2,0.5'), ('TU', '0.2,0.5'), ('JQ_AF_SK', '0.33,0.21')]), tag.format("C,T", "filter", "", format_dict, 0))
# 
#         format_dict = OrderedDict(zip("A:TAR:DP2".split(":"), "1:13,32:53".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('TAR', '13,32'), ('DP2', '53'), ('JQ_AF_SK', '0.6')]), tag.format("C", "filter", "", format_dict, 0))
# 
#         format_dict = OrderedDict(zip("A:TAR:DP2".split(":"), "1:13,50:53".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('TAR', '13,50'), ('DP2', '53'), ('JQ_AF_SK', '0.94')]), tag.format("C", "filter", "", format_dict, 0))
# 
#         format_dict = OrderedDict(zip("A:TAR:DP2".split(":"), "1:13,70:53".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('TAR', '13,70'), ('DP2', '53'), ('JQ_AF_SK', '1.00')]), tag.format("C", "filter", "", format_dict, 0))
# 
# class Strelka_DepthTagTestCase(unittest.TestCase):
#     def test_metaheader(self):
#         self.assertEqual('##FORMAT=<ID=JQ_DP_SK,Number=1,Type=Float,Description="Jacquard depth for Strelka (uses DP2 if available, otherwise uses ACGT tier2 depth),Source="Jacquard",Version=0.1>\n'.format(__version__), strelka.DepthTag().metaheader)
# 
#     def test_format_missingDP2AUTags(self):
#         tag = strelka.DepthTag()
#         format_param_string = "A:B"
#         format_value_string = "1:2"
#         format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("C", "filter", "", format_dict, 0))
# 
#     def test_format(self):
#         tag = strelka.DepthTag()
#         format_dict = OrderedDict(zip("A:DP2".split(":"), "1:42".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('DP2', '42'), ('JQ_DP_SK', '42')]), tag.format("C", "filter", "", format_dict, 0))
# 
#         format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:42,42:5,5:13,13:4,4".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('AU', '42,42'), ('CU', '5,5'), ('GU', '13,13'), ('TU', '4,4'), ('JQ_DP_SK', '64')]), tag.format("C", "filter", "", format_dict, 0))
# 
# class Strelka_SomaticTagTestCase(unittest.TestCase):
#     def test_metaheader(self):
#         self.assertEqual('##FORMAT=<ID=JQ_HC_SOM_SK,Number=1,Type=Integer,Description="Jacquard somatic status for Strelka: 0=non-somatic,1= somatic (based on PASS in FILTER column),Source="Jacquard",Version=0.1>\n'.format(__version__), strelka.SomaticTag().metaheader)
# 
#     def test_format_missingPASSTag(self):
#         tag = strelka.SomaticTag()
#         format_param_string = "A:B"
#         format_value_string = "1:2"
#         format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('B', '2'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "reject", "", format_dict, 0))
#         self.assertEqual(OrderedDict([('A', '1'), ('B', '2'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "PASS;foo", "", format_dict, 0))
# 
#     def test_format(self):
#         tag = strelka.SomaticTag()
#         format_dict = OrderedDict(zip("A:SS".split(":"), "1:2".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('SS', '2'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "PASS", "", format_dict, 0))
#         self.assertEqual(OrderedDict([('A', '1'), ('SS', '2'), ('JQ_HC_SOM_SK', '1')]), tag.format("C", "PASS", "", format_dict, 1))
# 
#         format_dict = OrderedDict(zip("A:SS".split(":"), "1:1".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('SS', '1'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "PASS", "", format_dict, 0))
#         self.assertEqual(OrderedDict([('A', '1'), ('SS', '1'), ('JQ_HC_SOM_SK', '1')]), tag.format("C", "PASS", "", format_dict, 1))

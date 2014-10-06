from collections import OrderedDict
import unittest

import jacquard.variant_callers.mutect as mutect
from jacquard.jacquard_utils import __version__, jq_dp_tag, jq_somatic_tag
from jacquard.vcf import VcfRecord 



class AlleleFreqTagTestCase(unittest.TestCase):

    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_MT,Number=A,Type=Float,Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA),Source="Jacquard",Version={0}>'.format(__version__), mutect._AlleleFreqTag().metaheader)
                
    def test_format_missingAFTag(self):
        tag = mutect._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        originalVcfRecord = VcfRecord(line)
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())
        
    def test_format_presentAFTag(self):
        tag = mutect._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3|0.567:SA.2:SA.3|0.834:SB.2:SB.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3:JQ_AF_MT|0.567:SA.2:SA.3:0.57|0.834:SB.2:SB.3:0.83\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
    def test_format_multAlt(self):
        tag = mutect._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3|0.5,0.8:SA.2:SA.3|0.7,0.6:SB.2:SB.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3:JQ_AF_MT|0.5,0.8:SA.2:SA.3:0.5,0.8|0.7,0.6:SB.2:SB.3:0.7,0.6\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
class DepthTagTestCase(unittest.TestCase):
    
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}MT,Number=1,Type=Float,Description="Jacquard depth for MuTect (based on DP),Source="Jacquard",Version={1}>'.format(jq_dp_tag, __version__), mutect._DepthTag().metaheader)
        
    def test_format_missingDPTag(self):
        tag = mutect._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        originalVcfRecord = VcfRecord(line)
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())
        
    def test_format_presentDPTag(self):
        tag = mutect._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3|2:SA.2:SA.3|4:SB.2:SB.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3:JQ_DP_MT|2:SA.2:SA.3:2|4:SB.2:SB.3:4\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
class SomaticTagTestCase(unittest.TestCase):
 
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1=somatic (based on SS FORMAT tag),Source="Jacquard",Version={1}>'.format(jq_somatic_tag, __version__), mutect._SomaticTag().metaheader)
        
    def test_format_missingSSTag(self):
        tag = mutect._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:"+jq_somatic_tag+"MT|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
    def test_format_presentSSTag(self):
        tag = mutect._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|SS:F2:F3|2:SA.2:SA.3|5:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|SS:F2:F3:"+jq_somatic_tag+"MT|2:SA.2:SA.3:1|5:SB.2:SB.3:0\n").replace('|',"\t")
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
        
class Mutect_TestCase(unittest.TestCase):
    
    def setUp(self):
        self.caller = mutect.Mutect()
        
    def test_validateInputFile_isValid(self):
        metaheaders = ["##MuTect"]
        self.assertTrue(self.caller.validate_input_file(metaheaders, "#column_header"))
    
    def test_validateInputFile_isNotValid(self):
        metaheaders = ["Foo"]
        self.assertFalse(self.caller.validate_input_file(metaheaders, "#column_header"))
    
    def test_addTags(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        self.caller.tags=[MockTag("mockTag", 42)]
        actual_line = self.caller.add_tags(VcfRecord(input_line))

        expected_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:mockTag|SA.1:SA.2:SA.3:42|SB.1:SB.2:SB.3:42\n".replace('|',"\t")

        self.assertEquals(expected_line, actual_line)
        
    def test_updateMetaheader(self):
        self.caller.tags=[MockTag("mockTag", 42, "##my_metaheader\n")]
        actual_metaheader = self.caller.get_new_metaheaders()

        self.assertEquals(["##my_metaheader\n"], actual_metaheader)


#     #TODO: (cgates/kmeng) Remove when switched to new VcfRecord
#     def test_format_rounds(self):
#         tag = mutect.AlleleFreqTag()
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.2".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.2'), ('JQ_AF_MT', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.20".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.20'), ('JQ_AF_MT', '0.20')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.204".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.204'), ('JQ_AF_MT', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.205".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.205'), ('JQ_AF_MT', '0.21')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.206".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.206'), ('JQ_AF_MT', '0.21')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:1.0".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '1.0'), ('JQ_AF_MT', '1.0')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:1.00".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '1.00'), ('JQ_AF_MT', '1.00')]), tag.format("alt", "filter", "", format_dict, 0))
# 
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.204,0.3807".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.204,0.3807'), ('JQ_AF_MT', '0.2,0.38')]), tag.format("alt", "filter", "", format_dict, 0))
#  
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.204,0.3807,0.2784".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.204,0.3807,0.2784'), ('JQ_AF_MT', '0.2,0.38,0.28')]), tag.format("alt", "filter", "", format_dict, 0))
# 
# class Mutect_DepthTagTestCase(unittest.TestCase):
#     def test_metaheader(self):
#         self.assertEqual('##FORMAT=<ID=JQ_DP_MT,Number=1,Type=Float,Description="Jacquard depth for MuTect (based on DP),Source="Jacquard",Version={0}>\n'.format(__version__), mutect.DepthTag().metaheader)
#                 
#     def test_format_missingDPTag(self):
#         tag = mutect.DepthTag()
#         format_param_string = "A:B"
#         format_value_string = "1:2"
#         format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("alt", "filter", "", format_dict, 0))
#                 
#     def test_format(self):
#         tag = mutect.DepthTag()
#         format_dict = OrderedDict(zip("A:DP".split(":"), "1:42".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('DP', '42'), ('JQ_DP_MT', '42')]), tag.format("alt", "filter", "", format_dict, 0))
# 
# class Mutect_SomaticTagTestCase(unittest.TestCase):
#     def test_metaheader(self):
#         self.assertEqual('##FORMAT=<ID=JQ_HC_SOM_MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag),Source="Jacquard",Version={0}>\n'.format(__version__), mutect.SomaticTag().metaheader)
#                 
#     def test_format_missingSSTag(self):
#         tag = mutect.SomaticTag()
#         format_param_string = "A:B"
#         format_value_string = "1:2"
#         format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('B', '2'), ("JQ_HC_SOM_MT", "0")]), tag.format("alt", "filter", "", format_dict, 0))
#                 
#     def test_format(self):
#         tag = mutect.SomaticTag()
#         format_dict = OrderedDict(zip( "A:SS".split(":"), "1:2".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('SS', '2'), ('JQ_HC_SOM_MT', '1')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip( "A:SS".split(":"), "1:1".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('SS', '1'), ('JQ_HC_SOM_MT', '0')]), tag.format("alt", "filter", "", format_dict, 0))

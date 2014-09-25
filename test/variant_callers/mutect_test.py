from collections import OrderedDict
import unittest

import jacquard.variant_callers.mutect as mutect
from jacquard.jacquard_utils import __version__
from jacquard.vcf_record import VcfRecord 


class Mutect_AlleleFreqTagTestCase(unittest.TestCase):
    def __init__(self):
        self.vcf_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_MT,Number=A,Type=Float,Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA),Source="Jacquard",Version={0}>\n'.format(__version__), mutect.AlleleFreqTag().metaheader)
                
                
    def test_format_missingAFTag(self):
        tag = mutect.AlleleFreqTag()
        vcfRecord = VcfRecord(self.vcf_line)
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format(vcfRecord))

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

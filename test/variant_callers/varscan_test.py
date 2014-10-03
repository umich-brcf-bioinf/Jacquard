# pylint: disable=C0103,C0301,R0904
from collections import OrderedDict,defaultdict
import os
import unittest

from jacquard.variant_callers import varscan
from jacquard.jacquard_utils import __version__

class VarScanTestCase(unittest.TestCase):
    def test_validateInputFile_valid(self):
        caller = varscan.Varscan() 
        input_file = ["##source=VarScan2", "foo"]
        name, valid = caller.validate_input_file(input_file)
        
        self.assertEquals("VarScan", name)
        self.assertEquals(1, valid)
        
    def test_validateInputFile_invalid(self):
        caller = varscan.Varscan()
        input_file = ["##bar", "foo"]
        name, valid = caller.validate_input_file(input_file)
         
        self.assertEquals("VarScan", name)
        self.assertEquals(0, valid)
#         
    def test_finalSteps(self):
        caller = varscan.Varscan()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        hc_candidates  = {os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_merged.Somatic.hc'): [
            os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_indel.Somatic.hc'), os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_snp.Somatic.hc')], os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_merged.LOH.hc'): [os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_indel.LOH.hc'), os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_snp.LOH.hc')], os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_merged.Germline.hc'): [os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_indel.Germline.hc'), os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_snp.Germline.hc')]}
        merge_candidates = {os.path.join(script_dir, '../reference_files/normalize_varscan_test/tiny_merged.vcf'): [os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_indel.vcf'), os.path.join(script_dir, '../reference_files/normalize_varscan_test/input/tiny_snp.vcf')]}
        output_dir = os.path.join(script_dir, "../reference_files/normalize_varscan_test/output")
        marked_as_hc = caller.final_steps(hc_candidates, merge_candidates, output_dir)
         
        self.assertEquals(['chr1^14397^CTGT^C'], marked_as_hc)
         
    def test_handleHCFiles(self):
        caller = varscan.Varscan()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        in_file = os.path.join(script_dir, "normalize_varscan_test", "input", "tiny_indel.Germline.hc")
        out_dir = os.path.join(script_dir, "normalize_varscan_test")
        hc_candidates = defaultdict(list)
        hc_candidates = caller.handle_hc_files(in_file, out_dir, hc_candidates)
         
        expected_hc_candidates = {os.path.join(script_dir, "normalize_varscan_test", "tiny_merged.Germline.hc"): [os.path.join(script_dir, "normalize_varscan_test", "input", "tiny_indel.Germline.hc")]}
        self.assertEquals(expected_hc_candidates, hc_candidates)
             
    def test_validateFileSet(self):
        caller = varscan.Varscan()
        all_keys = ["foo_merged.vcf", "foo_merged.Somatic.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
        sample_files = caller.validate_file_set(all_keys)
         
        self.assertEqual({'foo': ['.vcf', '.Somatic.hc', '.Germline.hc', '.LOH.hc']}, sample_files)
         
    def test_validateFileSet_invalid(self):
        caller = varscan.Varscan()
        all_keys = ["foo_merged.vcf", "foo_merged.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
         
        with self.assertRaises(SystemExit) as cm:
            sample_files = caller.validate_file_set(all_keys)
        self.assertEqual(cm.exception.code, 1)
     
    def test_validateFileSet_extraFiles(self):
        caller = varscan.Varscan()
        all_keys = ["foo_merged.bar.hc", "foo_merged.vcf", "foo_merged.Somatic.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
        sample_files = caller.validate_file_set(all_keys)
         
        self.assertEqual({'foo': ['.bar.hc', '.vcf', '.Somatic.hc', '.Germline.hc', '.LOH.hc']}, sample_files)
        
    def test_identifyHcVariants_VarScan(self):
        caller = varscan.Varscan()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/../reference_files/normalize_varscan_test/input/"
        hc_candidates = {"tiny_merged.vcf" : [input_dir + "tiny_indel.Germline.hc", input_dir + "tiny_indel.LOH.hc", input_dir + "tiny_indel.Somatic.hc", input_dir + "tiny_snp.Germline.hc", input_dir + "tiny_snp.LOH.hc", input_dir + "tiny_snp.Somatic.hc"]}
          
        hc_variants = caller.identify_hc_variants(hc_candidates)
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
          
    def test_markHcCandidates_VarScan(self):
        caller = varscan.Varscan()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/../reference_files/normalize_varscan_test/merged_vcf/"
        output_dir = script_dir + "/../reference_files/normalize_varscan_test/output/"
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
        merge_candidates = {input_dir + "tiny_merged.vcf" : [input_dir + "tiny_snp.vcf", input_dir + "tiny_indel.vcf"]}
          
        marked_as_hc = caller.mark_hc_variants(hc_variants, merge_candidates, output_dir)
          
        self.assertEqual(['chr1^14397^CTGT^C'], marked_as_hc)
         
class Varscan_AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_VS,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ)",Source="Jacquard",Version={0}>\n'.format(__version__), varscan.AlleleFreqTag().metaheader)
                 
    def test_format_missingAFTag(self):
        tag = varscan.AlleleFreqTag()
        info_string = ""
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("alt", "filter", info_string, format_dict, 0))
                 
    def test_format_rounds(self):
        tag = varscan.AlleleFreqTag()
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20%'), ('JQ_AF_VS', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
         
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.0%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.0%'), ('JQ_AF_VS', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
         
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.4%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.4%'), ('JQ_AF_VS', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
         
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.5%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.5%'), ('JQ_AF_VS', '0.21')]), tag.format("alt", "filter", "", format_dict, 0))
         
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.6%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.6%'), ('JQ_AF_VS', '0.21')]), tag.format("alt", "filter", "", format_dict, 0))
         
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:100%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '100%'), ('JQ_AF_VS', '1.0')]), tag.format("alt", "filter", "", format_dict, 0))
         
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:100.0%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '100.0%'), ('JQ_AF_VS', '1.0')]), tag.format("alt", "filter", "", format_dict, 0))
         
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.4%,38.075%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.4%,38.075%'), ('JQ_AF_VS', '0.2,0.38')]), tag.format("alt", "filter", "", format_dict, 0))
         
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.4%,38.075%,27.843%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.4%,38.075%,27.843%'), ('JQ_AF_VS', '0.2,0.38,0.28')]), tag.format("alt", "filter", "", format_dict, 0))
         
class Varscan_DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_DP_VS,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP)",Source="Jacquard",Version={0}>\n'.format(__version__), varscan.DepthTag().metaheader)
                 
    def test_format_missingDPTag(self):
        tag = varscan.DepthTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("alt", "filter", "", format_dict, 0))
                 
    def test_format(self):
        tag = varscan.DepthTag()
        format_dict = OrderedDict(zip("A:DP".split(":"), "1:42".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('DP', '42'), ('JQ_DP_VS', '42')]), tag.format("alt", "filter", "", format_dict, 0)) 


class Varscan_SomaticTagTestCase(unittest.TestCase):
    def setUp(self):
        self.tag = varscan.SomaticTag()

    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_HC_SOM_VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1= somatic (based on variant is high confidence, SOMATIC info field and if sample is TUMOR)",Source="Jacquard",Version={0}>\n'.format(__version__), varscan.SomaticTag().metaheader)

    def test_format_missingSSInfoTag(self):
        format_dict = OrderedDict([("A", "1"), ("B", "2")])
        actual_format_dict = self.tag.format("alt", "filter", "INFO",
                                             format_dict, 0)
        self.assertEqual(["A", "B", "JQ_HC_SOM_VS"], actual_format_dict.keys())
        self.assertEqual(["1", "2", "0"], actual_format_dict.values())

    def test_format_somaticTag0WhenSampleIsNormal(self):
        format_dict = OrderedDict([("A", "1")])
        sample_count = 0

        actual_format_dict = self.tag.format("alt", "filter", "INFO",
                                             format_dict, sample_count)
        expected_format_dict = OrderedDict([("A", "1"), ("JQ_HC_SOM_VS", "0")])
        self.assertEqual(expected_format_dict, actual_format_dict)

        actual_format_dict = self.tag.format("alt", "filter", "INFO;SS=2",
                                             format_dict, sample_count)
        expected_format_dict = OrderedDict([("A", "1"), ("JQ_HC_SOM_VS", "0")])
        self.assertEqual(expected_format_dict, actual_format_dict)


    def test_format_somaticTag1WhenSampleIsTumorAndSSAndHighConfidence(self):
        input_format_dict = {}
        sample_is_tumor = 1

        actual_format_dict = self.tag.format("alt", "filter",
                                             "INFO",
                                             input_format_dict, sample_is_tumor)
        expected_format_dict = OrderedDict([("JQ_HC_SOM_VS", "0")])
        self.assertEqual(expected_format_dict, actual_format_dict)

        actual_format_dict = self.tag.format("alt", "filter",
                                             "INFO;SS=2",
                                             input_format_dict, sample_is_tumor)
        expected_format_dict = OrderedDict([("JQ_HC_SOM_VS", "0")])
        self.assertEqual(expected_format_dict, actual_format_dict)

        actual_format_dict = self.tag.format("alt", "filter",
                                             "INFO;JQ_HC_VS",
                                             input_format_dict, sample_is_tumor)
        expected_format_dict = OrderedDict([("JQ_HC_SOM_VS", "0")])
        self.assertEqual(expected_format_dict, actual_format_dict)

        actual_format_dict = self.tag.format("alt", "filter",
                                             "INFO;JQ_HC_VS;SS=2",
                                             input_format_dict, sample_is_tumor)
        expected_format_dict = OrderedDict([("JQ_HC_SOM_VS", "1")])
        self.assertEqual(expected_format_dict, actual_format_dict)


    def test_format_rerunningUpdatesInsteadOfAdds(self):
        format_dict = OrderedDict([("JQ_HC_SOM_VS", "0")])
        sample_count = 1

        actual_format_dict = self.tag.format("alt", "filter",
                                             "INFO;SS=2;JQ_HC_VS",
                                             format_dict, sample_count)
        expected_format_dict = OrderedDict([("JQ_HC_SOM_VS", "1")])
        self.assertEqual(expected_format_dict, actual_format_dict)


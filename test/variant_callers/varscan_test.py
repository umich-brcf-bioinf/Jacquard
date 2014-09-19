from collections import OrderedDict,defaultdict
import glob
import os
from os import listdir
from StringIO import StringIO
import sys
from stat import *
import testfixtures
from testfixtures import TempDirectory
import unittest

import jacquard.variant_callers.varscan as varscan
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
         
class Varscan_AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_VS,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ),Source="Jacquard",Version={0}>\n'.format(__version__), varscan.AlleleFreqTag().metaheader)
                 
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
        self.assertEqual('##FORMAT=<ID=JQ_DP_VS,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP),Source="Jacquard",Version={0}>\n'.format(__version__), varscan.DepthTag().metaheader)
                 
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
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_HC_SOM_VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1= somatic (based on SOMATIC info tag and if sample is TUMOR),Source="Jacquard",Version={0}>\n'.format(__version__), varscan.SomaticTag().metaheader)
                 
    def test_format_missingSSInfoTag(self):
        tag = varscan.SomaticTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(["A", "B", "JQ_HC_SOM_VS"], tag.format("alt", "filter", "INFO", format_dict, 0).keys())
        self.assertEqual(["1", "2", "0"], tag.format("alt", "filter", "INFO", format_dict, 0).values())
                 
    def test_format(self):
        tag = varscan.SomaticTag()
        format_dict = OrderedDict([("A", "1")])
        self.assertEqual(OrderedDict([("A", "1"), ("JQ_HC_SOM_VS", "0")]), tag.format("alt", "filter", "INFO;SS=2", format_dict, 0))
        self.assertEqual(OrderedDict([("A", "1"), ("JQ_HC_SOM_VS", "1")]), tag.format("alt", "filter", "INFO;SS=2", format_dict, 1))
        self.assertEqual(OrderedDict([("A", "1"), ("JQ_HC_SOM_VS", "1")]), tag.format("alt", "filter", "INFO;SS=2;JQ_HC_SOM_VS", format_dict, 1))
         
        format_dict = OrderedDict([("A", "1")])
        self.assertEqual(OrderedDict([("A", "1",), ("JQ_HC_SOM_VS", "0")]), tag.format("alt", "filter", "INFO", format_dict, 0))
        self.assertEqual(OrderedDict([("A", "1"), ("JQ_HC_SOM_VS", "0")]), tag.format("alt", "filter", "INFO", format_dict, 1))

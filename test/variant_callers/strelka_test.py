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

import strelka
from jacquard_utils import __version__

class StrelkaTestCase(unittest.TestCase):
    def test_validateInputFile_valid(self):
        caller = strelka.Strelka()
        input_file = ["##source=strelka", "foo"]
        name, valid = caller.validate_input_file(input_file)
        
        self.assertEquals("Strelka", name)
        self.assertEquals(1, valid)
        
    def test_validateInputFile_invalid(self):
        caller = strelka.Strelka()
        input_file = ["##bar", "foo"]
        name, valid = caller.validate_input_file(input_file)
        
        self.assertEquals("Strelka", name)
        self.assertEquals(0, valid)
        
    def test_finalSteps(self):
        caller = strelka.Strelka()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        hc_candidates  = {os.path.join(script_dir, 'normalize_varscan_test/input/tiny_merged.Somatic.hc'): [os.path.join(script_dir, 'normalize_varscan_test/input/tiny_indel.Somatic.hc'), os.path.join(script_dir, 'normalize_varscan_test/input/tiny_snp.Somatic.hc')], os.path.join(script_dir, 'normalize_varscan_test/input/tiny_merged.LOH.hc'): [os.path.join(script_dir, 'normalize_varscan_test/input/tiny_indel.LOH.hc'), os.path.join(script_dir, 'normalize_varscan_test/input/tiny_snp.LOH.hc')], os.path.join(script_dir, 'normalize_varscan_test/input/tiny_merged.Germline.hc'): [os.path.join(script_dir, 'normalize_varscan_test/input/tiny_indel.Germline.hc'), os.path.join(script_dir, 'normalize_varscan_test/input/tiny_snp.Germline.hc')]}
        merge_candidates = {os.path.join(script_dir, 'normalize_varscan_test/output/tiny_merged.vcf'): [os.path.join(script_dir, 'normalize_varscan_test/input/tiny_indel.vcf'), os.path.join(script_dir, 'normalize_varscan_test/input/tiny_snp.vcf')]}
        output_dir = os.path.join(script_dir, "normalize_varscan_test/output")
        actual_merge_candidates = caller.final_steps(hc_candidates, merge_candidates, output_dir)
        
        self.assertEquals(merge_candidates, actual_merge_candidates)
        
    def test_handleHCFiles(self):
        caller = strelka.Strelka()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        in_file = os.path.join(script_dir, "normalize_varscan_test\input\tiny_indel.Germline.hc")
        out_dir = os.path.join(script_dir, "normalize_varscan_test\output")
        hc_candidates = defaultdict(list)
        hc_candidates = caller.handle_hc_files(in_file, out_dir, hc_candidates)
        
        self.assertEquals(defaultdict(list), hc_candidates)
            
    def test_validateFileSet(self):
        caller = strelka.Strelka()
        all_keys = ["foo_merged.vcf", "foo_merged.Somatic.hc", "foo_merged.Germline.hc", "foo_merged.LOH.hc"]
        sample_files = caller.validate_file_set(all_keys)
        
        self.assertEqual(None, sample_files)
        
class Strelka_AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_SK,Number=A,Type=Float,Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth),Source="Jacquard",Version={0}>\n'.format(__version__), strelka.AlleleFreqTag().metaheader)
                
    def test_format_missingAFTag(self):
        tag = strelka.AlleleFreqTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("C", "filter", "", format_dict, 0))
                
    def test_format_rounds(self):
        tag = strelka.AlleleFreqTag()
        format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:0.2,0.5:0.2,0.5:0.2,0.5:0.2,0.5:".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('AU', '0.2,0.5'), ('CU', '0.2,0.5'), ('GU', '0.2,0.5'), ('TU', '0.2,0.5'), ('JQ_AF_SK', '0.25')]), tag.format("C", "filter", "", format_dict, 0))
#         
        format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:0.2,0.6:0.2,0.6:0.2,0.5:0.2,0.5:".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('AU', '0.2,0.6'), ('CU', '0.2,0.6'), ('GU', '0.2,0.5'), ('TU', '0.2,0.5'), ('JQ_AF_SK', '0.27')]), tag.format("C", "filter", "", format_dict, 0))
#         
        format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:0.2,0.6:0.2,0.8:0.2,0.5:0.2,0.5:".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('AU', '0.2,0.6'), ('CU', '0.2,0.8'), ('GU', '0.2,0.5'), ('TU', '0.2,0.5'), ('JQ_AF_SK', '0.33,0.21')]), tag.format("C,T", "filter", "", format_dict, 0))
#         
        format_dict = OrderedDict(zip("A:TAR:DP2".split(":"), "1:13,32:53".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('TAR', '13,32'), ('DP2', '53'), ('JQ_AF_SK', '0.6')]), tag.format("C", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:TAR:DP2".split(":"), "1:13,50:53".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('TAR', '13,50'), ('DP2', '53'), ('JQ_AF_SK', '0.94')]), tag.format("C", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:TAR:DP2".split(":"), "1:13,70:53".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('TAR', '13,70'), ('DP2', '53'), ('JQ_AF_SK', '1.00')]), tag.format("C", "filter", "", format_dict, 0))

class Strelka_DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_DP_SK,Number=1,Type=Float,Description="Jacquard depth for Strelka (based on DP2),Source="Jacquard",Version=0.1>\n'.format(__version__), strelka.DepthTag().metaheader)
                
    def test_format_missingDP2AUTags(self):
        tag = strelka.DepthTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("C", "filter", "", format_dict, 0))
                
    def test_format(self):
        tag = strelka.DepthTag()
        format_dict = OrderedDict(zip("A:DP2".split(":"), "1:42".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('DP2', '42'), ('JQ_DP_SK', '42')]), tag.format("C", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:42,42:5,5:13,13:4,4".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('AU', '42,42'), ('CU', '5,5'), ('GU', '13,13'), ('TU', '4,4'), ('JQ_DP_SK', '64')]), tag.format("C", "filter", "", format_dict, 0))

class Strelka_SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_HC_SOM_SK,Number=1,Type=Integer,Description="Jacquard somatic status for Strelka: 0=non-somatic,1= somatic (based on PASS in FILTER column),Source="Jacquard",Version=0.1>\n'.format(__version__), strelka.SomaticTag().metaheader)
                
    def test_format_missingPASSTag(self):
        tag = strelka.SomaticTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "reject", "", format_dict, 0))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "PASS;foo", "", format_dict, 0))
                
    def test_format(self):
        tag = strelka.SomaticTag()
        format_dict = OrderedDict(zip( "A:SS".split(":"), "1:2".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('SS', '2'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "PASS", "", format_dict, 0))
        self.assertEqual(OrderedDict([('A', '1'), ('SS', '2'), ('JQ_HC_SOM_SK', '1')]), tag.format("C", "PASS", "", format_dict, 1))
        
        format_dict = OrderedDict(zip( "A:SS".split(":"), "1:1".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('SS', '1'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "PASS", "", format_dict, 0))
        self.assertEqual(OrderedDict([('A', '1'), ('SS', '1'), ('JQ_HC_SOM_SK', '1')]), tag.format("C", "PASS", "", format_dict, 1))

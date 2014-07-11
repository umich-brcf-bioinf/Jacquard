#!/usr/bin/python2.7
import unittest
from bin.tag_mutect import AlleleFreqTag, DepthTag, SomaticTag, LineProcessor

class AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_MT,Number=1,Type=Float, Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA).">', AlleleFreqTag().metaheader)
                
    def test_format_missingAFTag(self):
        tag = AlleleFreqTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(format_param_string, format_value_string))
                
    def test_format_rounds(self):
        tag = AlleleFreqTag()
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.2:0.2"), tag.format("A:FA", "1:0.2"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.20:0.20"), tag.format("A:FA", "1:0.20"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.204:0.2"), tag.format("A:FA", "1:0.204"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.205:0.21"), tag.format("A:FA", "1:0.205"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.206:0.21"), tag.format("A:FA", "1:0.206"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:1.0:1.0"), tag.format("A:FA", "1:1.0"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:1.00:1.00"), tag.format("A:FA", "1:1.00"))

class DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_DP_MT,Number=1,Type=Float, Description="Jacquard depth for MuTect (based on DP).">', DepthTag().metaheader)
                
    def test_format_missingDPTag(self):
        tag = DepthTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(format_param_string, format_value_string))
                
    def test_format(self):
        tag = DepthTag()
        self.assertEqual(("A:DP:JQ_DP_MT", "1:42:42"), tag.format("A:DP", "1:42"))

class SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_SOM_MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag).">', SomaticTag().metaheader)
                
    def test_format_missingSSTag(self):
        tag = SomaticTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(format_param_string, format_value_string))
                
    def test_format(self):
        tag = SomaticTag()
        self.assertEqual(("A:SS:JQ_SOM_MT", "1:2:1"), tag.format("A:SS", "1:2"))
        self.assertEqual(("A:SS:JQ_SOM_MT", "1:1:0"), tag.format("A:SS", "1:1"))

class LineProcessorTestCase(unittest.TestCase):
    def test_process_line_singleSample(self):
        tag = MockLowerTag()    
        processor = LineProcessor([tag])
        input_line = "chr1|42|.|ref|alt|qual|filter|INFO|A:B:C|X:Y:Z".replace("|", "\t")

        actual_line = processor.add_tags(input_line)
        actual_format_params, actual_format_values = actual_line.split("\t")[8:10]
        
        self.assertEqual("a:b:c", actual_format_params)
        self.assertEqual("x:y:z", actual_format_values)
        

class MockLowerTag():
    def format(self, params, values):
        return (params.lower(), values.lower())
    
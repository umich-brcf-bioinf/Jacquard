#!/usr/bin/python2.7
import unittest
from bin.tag_mutect import AlleleFreqTag

class AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_MT,Number=1,Type=Float, Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA).">', AlleleFreqTag().metaheader)
        
    def test_format(self):
        tag = AlleleFreqTag()
        format_param_string = "A:FA"
        format_value_string = "1:0.2"
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.2:0.2"), tag.format(format_param_string, format_value_string))
        
    def Xtest_format_missingAFTag(self):
        tag = AlleleFreqTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(format_param_string, format_value_string))
        
        
    def Xtest_format_rounds(self):
        tag = AlleleFreqTag()
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.234:0.23"), tag.format("A:FA", "1:0.234"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.234:0.23"), tag.format("A:FA", "1:0.234"))

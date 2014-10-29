import os 
import unittest

import jacquard.utils as utils
from jacquard.expand2 import _parse_meta_headers, _append_format_tags_to_samples, _get_headers

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

class MockVcfReader(object):
    def __init__(self, input_filepath="vcfName", metaheaders=["##metaheaders"], column_header="#header"):
        self.input_filepath = input_filepath
        self.metaheaders = metaheaders
        self.column_header = column_header
        self.opened = False
        self.closed = False

    def open(self):
        self.opened = True

    def vcf_records(self):
        return iter(["foo"])
    
    def close(self):
        self.closed = True

class ExpandTestCase(unittest.TestCase):
    def test_parse_meta_headers(self):
        meta_headers = ['##ALT=<ID=DEL,Description="Deletion">', 
                        '##INFO=<ID=AC,Number=.,Type=Integer,Description="Alternate Allele Count">]',
                        '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele"', 
                        '##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">',
                        '##RUNTIME_ARG=allele freq. cutoff: 5']
        (info_fields, format_tags) = _parse_meta_headers(meta_headers)
        
        self.assertEquals(["AA", "AC"], info_fields)
        self.assertEquals(["SP"], format_tags)
        
    def test_parse_meta_headers_missing(self):
        meta_headers = ['##ALT=<ID=DEL,Description="Deletion">', 
                        '##INFO=<ID=AC,Number=.,Type=Integer,Description="Alternate Allele Count">]',
                        '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele"',
                        '##RUNTIME_ARG=allele freq. cutoff: 5']

        self.assertRaisesRegexp(utils.JQException,
                                    "Unable to parse meta_headers for INFO and/or " + 
                                    "FORMAT fields. Review input and try again.",
                                    _parse_meta_headers,
                                    meta_headers)
        
    def test_append_format_tags_to_samples(self):
        format_tags = ["bar", "foo"] 
        samples = ["sampleA", "sampleB"]
        actual = _append_format_tags_to_samples(format_tags, samples)
        
        expected = ["bar|sampleA", "foo|sampleA", "bar|sampleB", "foo|sampleB"]
        
        self.assertEquals(expected, actual)
        
    def test_get_headers(self):
        meta_headers = ['##ALT=<ID=DEL,Description="Deletion">', 
                        '##INFO=<ID=AC,Number=.,Type=Integer,Description="Alternate Allele Count">]',
                        '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele"', 
                        '##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">',
                        '##RUNTIME_ARG=allele freq. cutoff: 5']
        
        mock_reader = MockVcfReader(metaheaders=meta_headers, column_header="CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\tsampleB")
        actual = _get_headers(mock_reader)
        
        expected = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tFORMAT\tAA\tAC\tSP|sampleA\tSP|sampleB"
        
        self.assertEquals(expected, actual)
        
    
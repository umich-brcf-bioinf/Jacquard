from collections import OrderedDict
import os 
import unittest

import jacquard.utils as utils
from jacquard.expand2 import _parse_meta_headers, _append_format_tags_to_samples, _get_headers, _write_vcf_records, _disambiguate_column_names

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

class MockVcfReader(object):
    def __init__(self, input_filepath="vcfName", metaheaders=["##metaheaders"], column_header="#header", content = ["foo"]):
        self.input_filepath = input_filepath
        self.metaheaders = metaheaders
        self.column_header = column_header
        self.opened = False
        self.closed = False
        self.content = content
        
    def open(self):
        self.opened = True

    def vcf_records(self):
        for content in self.content:
            yield MockVcfRecord(content)
    
    def close(self):
        self.closed = True
        
class MockVcfRecord(object):
    def __init__(self, content):
        self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format = content[0:9]
        self.samples = content[9:]
        
        tags = self.format.split(":")
        self.sample_dict = {}
        for i, sample in enumerate(self.samples):
            values = sample.split(":")
            self.sample_dict[i] = OrderedDict(zip(tags, values))
            
    def get_info_dict(self):
        info_dict = {}
        
        for key_value in self.info.split(";"):
            if "=" in key_value:
                key,value = key_value.split("=")
                info_dict[key] = value
            else:
                info_dict[key_value] = key_value
        
        return info_dict
    
class MockFileWriter(object):
    def __init__(self):
        self.written = []
        
    def write(self, text):
        self.written.append(text)

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
        
    def test_disambiguate_column_names(self):
        actual = _disambiguate_column_names(["CHROM", "POS", "ID", "REF"], ["HOM", "AA", "DP", "GT"])
        expected = ["HOM", "AA", "DP", "GT"]
        
        self.assertEquals(expected, actual)
        
        actual = _disambiguate_column_names(["CHROM", "POS", "ID", "REF"], ["HOM", "AA", "ID", "DP", "GT"])
        expected = ["INFO_HOM", "INFO_AA", "INFO_ID", "INFO_DP", "INFO_GT"]
        
        self.assertEquals(expected, actual)
        
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
        
        expected = (["CHROM","POS","ID","REF","ALT","QUAL","FILTER"],["AA","AC"],["SP|sampleA","SP|sampleB"])
        
        self.assertEquals(expected, actual)
        
    def test_write_vcf_records(self):
        column_header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA"
        
        mock_vcf_reader = MockVcfReader(content=[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","tag1=val1;tag3=val3;tag4","FOO:BAR","42:1"]], column_header=column_header)
        mock_file_writer = MockFileWriter()
        
        info_header = ["tag1", "tag2", "tag3", "tag4"]
        format_sample_header = ["BAR|sampleA", "FOO|sampleA"]
        
        _write_vcf_records(mock_vcf_reader, mock_file_writer, info_header, format_sample_header)
        
        actual = mock_file_writer.written
        expected = ["CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tval1\t\tval3\ttag4\t1\t42\n"]
        
        self.assertEquals(expected, actual)
        

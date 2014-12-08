#pylint: disable=C0111,C0301
from argparse import Namespace
from collections import OrderedDict
import os
from testfixtures import TempDirectory
import unittest

import jacquard.merge2 as merge2

class MockVcfReader(object):
    def __init__(self,
                 input_filepath="vcfName",
                 metaheaders=["##metaheaders"],
                 column_header="#header",
                 content=["foo"]):
        self.file_name = input_filepath
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
        self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, \
            self.filter, self.info, self.format = content[0:9]
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
                key, value = key_value.split("=")
                info_dict[key] = value
            else:
                info_dict[key_value] = key_value

        return info_dict

class MockFileWriter(object):
    def __init__(self):
        self.written = []

    def write(self, text):
        self.written.append(text)

class MergeTestCase(unittest.TestCase):
    def test_get_metaheaders(self):
        meta_headers = ['##fileformat=VCFv4.2',
                        '##jacquard.version=0.21',
                        '##FORMAT=<ID=JQ_MT_AF,Number=A,Type=Float,Description="foo",Source="Jacquard",Version=0.21>',
                        '##contig=<ID=chr1,length=249350621,assembly=hg19']
        mock_vcf_reader = MockVcfReader(input_filepath = "P1.vcf",
                                        metaheaders = meta_headers,
                                        column_header = 'CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR')

        actual_meta_headers = merge2._produce_merged_metaheaders(mock_vcf_reader, [], 1)

        meta_headers.append("##jacquard.merge.file1=P1.vcf(['NORMAL', 'TUMOR'])")
        meta_headers.append('##INFO=<ID=JQ_MULT_ALT_LOCUS,Number=0,Type=Flag,Description="dbSNP Membership",Source="Jacquard",Version="0.21">')

        self.assertEquals(meta_headers, actual_meta_headers)

#TODO: fix
    def xtest_execute(self):
        vcf_content1 = '''##fileformat=VCFv4.2
##source=strelka
##jacquard.version=0.21
##FORMAT=<ID=JQ_MT_AF,Number=A,Type=Float,Description="foo",Source="Jacquard",Version=0.21>
##contig=<ID=chr1,length=249350621,assembly=hg19
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
'''
        vcf_content2 = '''##fileformat=VCFv4.2
##source=strelka
##jacquard.version=0.21
##FORMAT=<ID=JQ_MT_AF,Number=A,Type=Float,Description="foo",Source="Jacquard",Version=0.21>
##contig=<ID=chr1,length=249350621,assembly=hg19
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
'''
        vcf_content1 = vcf_content1.replace('|', "\t")
        vcf_content2 = vcf_content2.replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf", vcf_content1)
            input_dir.write("B.vcf", vcf_content2)

            args = Namespace(input=input_dir.path,
                             output=os.path.join(output_dir.path, "merged.vcf"))

            merge2.execute(args, [])

            output_dir.check("merged.vcf")
            file_content1 = output_dir.read("merged.vcf")

            self.assertTrue("##jacquard.merge.file1=A.vcf(['NORMAL', 'TUMOR'])" in file_content1)

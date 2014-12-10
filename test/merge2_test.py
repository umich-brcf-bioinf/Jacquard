#pylint: disable=C0111,C0301
from argparse import Namespace
from collections import OrderedDict
import glob
import os
from testfixtures import TempDirectory
import unittest

import jacquard.merge2 as merge2
import jacquard.utils as utils
import jacquard.vcf as vcf

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

    def test_functional_merge2(self):
        with TempDirectory() as output_dir:
            filedir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(filedir,
                                          "functional_tests",
                                          "04_merge2")
            input_dir = os.path.join(module_testdir,"input")
            args = Namespace(input=input_dir,
                             output=os.path.join(output_dir.path,"merged.vcf"),
                             allow_inconsistent_sample_sets=False)

            execution_context = ["##jacquard.version={0}".format(utils.__version__),
                                 "##jacquard.command=merge_test",
                                 "##jacquard.cwd=foo"]

            merge2.execute(args, execution_context)

            output_file = glob.glob(os.path.join(output_dir.path, "*.vcf"))[0]

            actual_file = vcf.FileReader(output_file)
            actual_file.open()
            actual = []
            for line in actual_file.read_lines():
                actual.append(line)
            actual_file.close()

            module_outdir = os.path.join(module_testdir,"benchmark")
            output_file = glob.glob(os.path.join(module_outdir, "*.vcf"))[0]

            expected_file = vcf.FileReader(os.path.join(module_outdir,output_file))
            expected_file.open()
            expected = []
            for line in expected_file.read_lines():
                expected.append(line)
            expected_file.close()
            print actual
            self.assertEquals(len(expected), len(actual))
            self.assertEquals(94, len(actual))

            for i in xrange(len(expected)):
                if expected[i].startswith("##jacquard.cwd="):
                    self.assertTrue(actual[i].startswith("##jacquard.cwd="))
                elif expected[i].startswith("##jacquard.command="):
                    self.assertTrue(actual[i].startswith("##jacquard.command="))
                else:
                    self.assertEquals(expected[i].rstrip(), actual[i].rstrip())


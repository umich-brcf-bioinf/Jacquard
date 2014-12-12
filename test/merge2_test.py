#pylint: disable=C0111,C0301,W0212
from argparse import Namespace
from collections import OrderedDict, defaultdict
import glob
import os
from testfixtures import TempDirectory
import unittest

import jacquard.merge2 as merge2
import jacquard.utils as utils
import jacquard.vcf as vcf
import test_case as test_case

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
            self.filter, self.info, self.format = content.split("\t")[0:9]
        self.samples = content.split("\t")[9:]

        tags = self.format.split(":")
        self.format_set = tags

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

    def asText(self):
        stringifier = [self.chrom, self.pos, self.id, self.ref, self.alt,
                   self.qual, self.filter, self.info,
                   ":".join(self.format_set)]

        for key in self.sample_dict:
            stringifier.append(":".join(self.sample_dict[key].values()))

        return "\t".join(stringifier) + "\n"

class MockFileWriter(object):
    def __init__(self):
        self.written = []

    def write(self, text):
        self.written.append(text)

class MergeTestCase(unittest.TestCase):
    def test_produce_merged_metaheaders(self):
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

    def test_add_to_coordinate_dict(self):
        mock_readers = [MockVcfReader(content=["chr1\t1\t.\tA\tC\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR",
                                               "chr2\t12\t.\tA\tG\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR"]),
                        MockVcfReader(content=["chr42\t16\t.\tG\tC\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR"])]
        coordinate_dict = defaultdict(list)
        for mock_reader in mock_readers:
            coordinate_dict = merge2._add_to_coordinate_dict(mock_reader, coordinate_dict)

        expected_coordinates = {"chr1^1": ["chr1^1^A^C"],
                                "chr2^12": ["chr2^12^A^G"],
                                "chr42^16": ["chr42^16^G^C"]}
        self.assertEquals(expected_coordinates, coordinate_dict)

    def test_add_to_coordinate_dict_multAlts(self):
        mock_readers = [MockVcfReader(content=["chr1\t1\t.\tA\tC\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR",
                                               "chr2\t12\t.\tA\tG\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR"]),
                        MockVcfReader(content=["chr1\t1\t.\tA\tT\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR",
                                               "chr42\t16\t.\tG\tC\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR"])]
        coordinate_dict = defaultdict(list)
        for mock_reader in mock_readers:
            coordinate_dict = merge2._add_to_coordinate_dict(mock_reader, coordinate_dict)

        expected_coordinates = {"chr1^1": ["chr1^1^A^C", "chr1^1^A^T"],
                                "chr2^12": ["chr2^12^A^G"],
                                "chr42^16": ["chr42^16^G^C"]}
        self.assertEquals(expected_coordinates, coordinate_dict)

    def test_sort_coordinate_dict(self):
        coordinate_dict = {"chr1^1": ["chr1^1^A^C", "chr1^1^A^T"],
                           "chr12^24": ["chr12^24^A^G"],
                           "chr4^16": ["chr4^16^G^C"]}

        sorted_coordinates = merge2._sort_coordinate_dict(coordinate_dict)
        expected_coordinates = OrderedDict([("chr1^1", ["chr1^1^A^C", "chr1^1^A^T"]),
                                            ("chr4^16", ["chr4^16^G^C"]),
                                            ("chr12^24", ["chr12^24^A^G"])])

        self.assertEquals(expected_coordinates, sorted_coordinates)

    def test_sort_coordinate_dict_multAlts(self):
        coordinate_dict = {"chr1^1": ["chr1^1^A^C"],
                           "chr12^24": ["chr12^24^A^G"],
                           "chr4^16": ["chr4^16^G^C"]}

        sorted_coordinates = merge2._sort_coordinate_dict(coordinate_dict)
        expected_coordinates = OrderedDict([("chr1^1", ["chr1^1^A^C"]),
                                            ("chr4^16", ["chr4^16^G^C"]),
                                            ("chr12^24", ["chr12^24^A^G"])])

        self.assertEquals(expected_coordinates, sorted_coordinates)

    def test_get_record_sample_data(self):
        mock_vcf_reader = MockVcfReader(content=["chr1\t1\t.\tA\tC\t.\t.\tINFO\tDP:AF\t42:0.23\t25:0.77",
                                                 "chr2\t12\t.\tA\tG\t.\t.\tINFO\tDP:AF\t41:0.67\t56:0.21"])
        format_tags = ["DP", "AF", "foo"]

        records = []
        for mock_vcf_record in mock_vcf_reader.vcf_records():
            records.append(mock_vcf_record)
        samples = merge2._get_record_sample_data(records[0], format_tags)

        expected_samples = {0: OrderedDict([("DP", "42"), ("AF", "0.23")]),
                            1: OrderedDict([("DP", "25"), ("AF", "0.77")])}
        self.assertEquals(expected_samples, samples)

    def test_write_variants(self):
        mock_reader = MockVcfReader(content=["chr1\t1\t.\tA\tC\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR",
                                             "chr2\t12\t.\tA\tG\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR"])
        mock_writer = MockFileWriter()
        merge2._write_variants(mock_reader, mock_writer, ["FORMAT"],["chr2^12^A^G"])

        expected = ["chr2\t12\t.\tA\tG\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR\n"]
        self.assertEquals(expected, mock_writer.written)

    def test_write_variants_multAlts(self):
        mock_reader = MockVcfReader(content=["chr1\t1\t.\tA\tC\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR",
                                             "chr1\t1\t.\tA\tT\t.\t.\t.\tFORMAT\tNORMAL\tTUMOR",
                                             "chr2\t12\t.\tA\tG\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR"])
        mock_writer = MockFileWriter()
        merge2._write_variants(mock_reader, mock_writer, ["FORMAT"], ["chr1^1^A^C", "chr1^1^A^T"])

        expected = ["chr1\t1\t.\tA\tC\t.\t.\tINFO;JQ_MULT_ALT_LOCUS\tFORMAT\tNORMAL\tTUMOR\n",
                    "chr1\t1\t.\tA\tT\t.\t.\tJQ_MULT_ALT_LOCUS\tFORMAT\tNORMAL\tTUMOR\n"]

        self.assertEquals(expected, mock_writer.written)

    def test_write_variants_missingCoordinate(self):
        mock_reader = MockVcfReader(content=["chr1\t1\t.\tA\tC\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR",
                                             "chr2\t12\t.\tA\tG\t.\t.\tINFO\tFORMAT\tNORMAL\tTUMOR"])
        mock_writer = MockFileWriter()
        merge2._write_variants(mock_reader, mock_writer, ["FORMAT"],["chr10^12^A^G"])

        expected = []
        self.assertEquals(expected, mock_writer.written)

class Merge2FunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_merge2(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "04_merge2")
            input_dir = os.path.join(module_testdir, "input")
            output_file = os.path.join(output_dir.path, "merged.vcf")

            command = ["merge2", input_dir, output_file, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

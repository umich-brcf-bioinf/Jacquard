#pylint: disable=missing-docstring,line-too-long,too-many-public-methods
#pylint: disable=too-few-public-methods,too-many-instance-attributes
#pylint: disable=too-many-arguments,invalid-name
from collections import OrderedDict
import os
from sets import Set
from testfixtures import TempDirectory
import unittest

import jacquard.merge2 as merge2
import jacquard.vcf as vcf
import test_case as test_case
from jacquard.vcf import VcfRecord

class MockVcfReader(object):
    def __init__(self,
                 input_filepath="vcfName",
                 metaheaders=["##metaheaders"],
                 column_header="#header",
                 content=["foo"],
                 records=None,
                 samples=None):
        self.content = content
        if records is None:
            self.records = [MockVcfRecord.parse_record(line) for line in self.content]
        else:
            self.records = records
        if samples is None:
            self.samples = []
        else:
            self.samples = samples
        self.file_name = input_filepath
        self.input_filepath = input_filepath
        self.metaheaders = metaheaders
        self.column_header = column_header
        self.opened = False
        self.closed = False

    def open(self):
        self.opened = True

    def vcf_records(self):
        for record in self.records:
            yield record

    def close(self):
        self.closed = True

class MockVcfRecord(object):
    @classmethod
    def parse_record(cls, vcf_line):
        vcf_fields = vcf_line.rstrip().split("\t")
        chrom, pos, rid, ref, alt, qual, rfilter, info, rformat \
                = vcf_fields[0:9]
        samples = vcf_fields[9:]
        return MockVcfRecord(chrom, pos, ref, alt, rid, qual, rfilter, info,
                             rformat, samples)

    def __init__(self, chrom, pos, ref, alt,
                 vcf_id=".", qual=".", vcf_filter=".", info=".", vcf_format=".",
                 samples=None):
        self.chrom = chrom
        self.pos = pos
        self.id = vcf_id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = vcf_filter
        self.info = info
        self.format = vcf_format
        if samples is None:
            self.samples = []
        else:
            self.samples = samples

        tags = self.format.split(":")
        self.format_set = tags

        self.sample_dict = {}
        for i, sample in enumerate(self.samples):
            values = sample.split(":")
            self.sample_dict[i] = OrderedDict(zip(tags, values))

    def get_empty_record(self):
        return MockVcfRecord(self.chrom, self.pos, self.ref, self.alt)

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

    def __eq__(self, other):
        return ("^".join([self.chrom,
                          self.pos,
                          self.ref,
                          self.alt]) == other)

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
        mock_vcf_reader = MockVcfReader(input_filepath="P1.vcf",
                                        metaheaders=meta_headers,
                                        column_header='CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR',
                                        records=[])

        actual_meta_headers = merge2._produce_merged_metaheaders(mock_vcf_reader, [], 1)

        meta_headers.append("##jacquard.merge.file1=P1.vcf(['NORMAL', 'TUMOR'])")
        meta_headers.append('##INFO=<ID=JQ_MULT_ALT_LOCUS,Number=0,Type=Flag,Description="dbSNP Membership",Source="Jacquard",Version="0.21">')

        self.assertEquals(meta_headers, actual_meta_headers)

    def test_extract_format_ids(self):
        meta_headers = ['##fileformat=VCFv4.2',
                        '##FORMAT=<Number=A,Type=Float,ID=JQ_MT_AF,Description="foo",Source="Jacquard",Version=0.21>',
                        '##FORMAT=<ID=JQ_VS_AF,Number=A,Type=Float,Description="foo",Source="Jacquard",Version=0.21>',
                        '##FORMAT=<Number=A,Type=Float,Description="foo",Source="Jacquard",Version=0.21,ID=JQ_SK_AF>',
                        '##contig=<ID=chr1,length=249350621,assembly=hg19']

        actual_format_tags = merge2._extract_format_ids(meta_headers)
        expected_format_tags = Set(["JQ_MT_AF", "JQ_VS_AF", "JQ_SK_AF"])

        self.assertEquals(expected_format_tags, actual_format_tags)

    def test_extract_format_ids_malformedMetaHeaders(self):
        meta_headers = ['##FORMAT=<Number=A>',
                        '##FORMAT=<ID1=JQ_MT_AF>',
                        '##FORMAT=<ID=JQ_SK1_AF,ID=JQ_SK2_AF>',
                        '##FORMAT=<ID="JQ_SK_AF">',
                        '##FORMAT=<ID=JQ_VS_AF>',
                        '##FORMAT=<ID=JQ_VS_AF>']

        actual_format_tags = merge2._extract_format_ids(meta_headers)
        expected_format_tags = Set(["JQ_SK2_AF", "JQ_VS_AF"])

        self.assertEquals(expected_format_tags, actual_format_tags)


    def test_build_coordinates(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "G", "id=2")
        fileBrec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileArec1, fileArec2]),
                        MockVcfReader(records=[fileBrec1, fileBrec2])]

        actual = merge2._build_coordinates(mock_readers)

        expected = [fileArec1, fileArec2, fileBrec2]
        self.assertEquals(expected, actual)

    def test_build_coordinates_multAltsEmpty(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "G", "id=2")
        fileBrec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileArec1, fileArec2]),
                        MockVcfReader(records=[fileBrec1, fileBrec2])]

        actual = merge2._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = []
        self.assertEquals(expected, actual_multalts)

    def test_build_coordinates_multAlts(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "T", "id=2")
        fileBrec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileArec1, fileArec2]),
                        MockVcfReader(records=[fileBrec1, fileBrec2])]

        actual = merge2._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = [fileArec2, fileBrec1]
        self.assertEquals(expected, actual_multalts)

    def Xtest_get_sample_tag_values(self):
        OD = OrderedDict
        samples1 = OD({"SA": OD({"foo":"A1", "bar":"A2"}),
                       "SB": OD({"foo":"B1", "bar":"B2"})})
        samples2 = OD({"SC": OD({"foo":"C1", "bar":"C2"}),
                       "SD": OD({"foo":"D1", "bar":"D2"})})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)
        record3 = None
        reader1 = MockBufferedReader([record1])
        reader2 = MockBufferedReader([record2])
        reader3 = MockBufferedReader([record3])
        buffered_readers = [reader1, reader2, reader3]
        merged_record = VcfRecord("chr1", "1", "A", "C")

        actual_sample_tag_values = merge2._get_sample_tag_values(buffered_readers, merged_record)

        self.assertEqual(set(["SA", "SB", "SC", "SD"]), set(actual_sample_tag_values.keys()))

    def Xtest_sort_coordinate_set_multAlts(self):
        coordinate_dict = {"chr1^1^A": ["chr1^1^A^C"],
                           "chr12^24^A": ["chr12^24^A^G"],
                           "chr4^16^G": ["chr4^16^G^C"]}

        sorted_coordinates = merge2._sort_coordinate_set(coordinate_dict)
        expected_coordinates = OrderedDict([("chr1^1^A", ["chr1^1^A^C"]),
                                            ("chr4^16^G", ["chr4^16^G^C"]),
                                            ("chr12^24^A", ["chr12^24^A^G"])])

        self.assertEquals(expected_coordinates, sorted_coordinates)

    def test_get_record_sample_data(self):
        mock_vcf_reader = MockVcfReader(content=["chr1\t1\t.\tA\tC\t.\t.\tINFO\tDP:AF\t42:0.23\t25:0.77",
                                                 "chr2\t12\t.\tA\tG\t.\t.\tINFO\tDP:AF\t41:0.67\t56:0.21"])
        format_tags = ["DP", "AF", "foo"]

        records = []
        for mock_vcf_record in mock_vcf_reader.vcf_records():
            records.append(mock_vcf_record)
        samples = merge2._get_record_sample_data(records[0], format_tags)

        expected_samples = {0: OrderedDict([("DP", "42"), ("AF", "0.23"), ("foo", ".")]),
                            1: OrderedDict([("DP", "25"), ("AF", "0.77"), ("foo", ".")])}
        self.assertEquals(expected_samples, samples)

#pylint: disable=unused-argument
class MockBufferedReader(object):
    def __init__(self, vcf_records):
        self.vcf_records_iter = iter(vcf_records)

    def get_if_equals(self, requested_record):
        return self.vcf_records_iter.next()

class Merge2FunctionalTestCase(test_case.JacquardBaseTestCase):
    def xtest_merge2(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "04_merge2")
            input_dir = os.path.join(module_testdir, "input")
            output_file = os.path.join(output_dir.path, "merged.vcf")

            command = ["merge2", input_dir, output_file, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)


class BufferedReaderTestCase(test_case.JacquardBaseTestCase):
    def test_get_sample_info_emptyIfNotRequestedCoordinateDoesNotMatch(self):
        rec1 = VcfRecord.parse_record("chr2\t2\t.\tA\tG\t.\tPASS\tINFO\tDP\t42\t16", ["SA", "SB"])

        mock_reader = MockVcfReader(input_filepath="fileA.vcf",
                                    records=[rec1],
                                    samples=["NORMAL", "TUMOR"])
        buffered_reader = merge2.BufferedReader(mock_reader)

        input_coordinate = VcfRecord("chr1", "1", "X", "X")
        self.assertEquals({}, buffered_reader.get_sample_info(input_coordinate))

    def Xtest_get_sample_info(self):
        rec1 = VcfRecord.parse_record("chr2\t2\t.\tA\tG\t.\tPASS\tINFO\tDP\t42\t16", ["SA", "SB"])

        mock_reader = MockVcfReader(input_filepath="fileA.vcf",
                                    records=[rec1],
                                    samples=["NORMAL", "TUMOR"])
        buffered_reader = merge2.BufferedReader(mock_reader)

        input_coordinate = VcfRecord("chr2", "2", "A", "G")
        actual_sample_info = buffered_reader.get_sample_info(input_coordinate)

        self.assertEquals(["fileA.vcf|NORMAL", "fileA.vcf|TUMOR"], sorted(actual_sample_info.keys()))
        self.assertEquals(OrderedDict({"DP": "42"}), actual_sample_info["fileA.vcf|NORMAL"])
        self.assertEquals(OrderedDict({"DP": "16"}), actual_sample_info["fileA.vcf|TUMOR"])

class GenericBufferedReaderTestCase(test_case.JacquardBaseTestCase):
    def test_get_sample_info_advancesCurrentElementWhenMatched(self):
        reader = [1, 5, 10, 15]
        buffered_reader = merge2.GenericBufferedReader(iter(reader))

        self.assertEquals(None, buffered_reader.get_if_equals(0))
        self.assertEquals(None, buffered_reader.get_if_equals(5))
        self.assertEquals(1, buffered_reader.get_if_equals(1))
        self.assertEquals(None, buffered_reader.get_if_equals(1))
        self.assertEquals(5, buffered_reader.get_if_equals(5))
        self.assertEquals(10, buffered_reader.get_if_equals(10))
        self.assertEquals(15, buffered_reader.get_if_equals(15))
        self.assertEquals(None, buffered_reader.get_if_equals(15))
        self.assertEquals(None, buffered_reader.get_if_equals(42))
        
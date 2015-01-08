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
                 sample_names=None):
        self.content = content
        if records is None:
            self.records = [MockVcfRecord.parse_record(line) for line in self.content]
        else:
            self.records = records
        if sample_names is None:
            self.sample_names = []
        else:
            self.sample_names = sample_names
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

#pylint: disable=unused-argument
class MockBufferedReader(object):
    def __init__(self, vcf_records):
        self.vcf_records_iter = iter(vcf_records)

    def get_if_equals(self, requested_record):
        return self.vcf_records_iter.next()
    
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

#     def test_extract_format_ids(self):
#         meta_headers = ['##fileformat=VCFv4.2',
#                         '##FORMAT=<Number=A,Type=Float,ID=JQ_MT_AF,Description="foo",Source="Jacquard",Version=0.21>',
#                         '##FORMAT=<ID=JQ_VS_AF,Number=A,Type=Float,Description="foo",Source="Jacquard",Version=0.21>',
#                         '##FORMAT=<Number=A,Type=Float,Description="foo",Source="Jacquard",Version=0.21,ID=JQ_SK_AF>',
#                         '##contig=<ID=chr1,length=249350621,assembly=hg19']
# 
#         actual_format_tags = merge2._extract_format_ids(meta_headers)
#         expected_format_tags = Set(["JQ_MT_AF", "JQ_VS_AF", "JQ_SK_AF"])
# 
#         self.assertEquals(expected_format_tags, actual_format_tags)

#     def test_extract_format_ids_malformedMetaHeaders(self):
#         meta_headers = ['##FORMAT=<Number=A>',
#                         '##FORMAT=<ID1=JQ_MT_AF>',
#                         '##FORMAT=<ID=JQ_SK1_AF,ID=JQ_SK2_AF>',
#                         '##FORMAT=<ID="JQ_SK_AF">',
#                         '##FORMAT=<ID=JQ_VS_AF>',
#                         '##FORMAT=<ID=JQ_VS_AF>']
# 
#         actual_format_tags = merge2._extract_format_ids(meta_headers)
#         expected_format_tags = Set(["JQ_SK2_AF", "JQ_VS_AF"])
# 
#         self.assertEquals(expected_format_tags, actual_format_tags)


    def test_build_coordinates(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "G", "id=2")
        fileBrec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileArec1, fileArec2]),
                        MockVcfReader(records=[fileBrec1, fileBrec2])]

        actual_coordinates = merge2._build_coordinates(mock_readers)[0]

        expected = [fileArec1, fileArec2, fileBrec2]
        self.assertEquals(expected, actual_coordinates)

    def test_build_coordinates_sortsSampleNames(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "G")

        mock_readers = [MockVcfReader(records=[fileArec1], sample_names=["SB", "SD"]),
                        MockVcfReader(records=[fileBrec1], sample_names=["SA", "SC"])]

        actual_samples = merge2._build_coordinates(mock_readers)[1]

        self.assertEquals(["SA", "SB", "SC", "SD"], actual_samples)

    def test_build_coordinates_multAltsEmpty(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "G", "id=2")
        fileBrec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileArec1, fileArec2]),
                        MockVcfReader(records=[fileBrec1, fileBrec2])]

        actual_coordinates = merge2._build_coordinates(mock_readers)[0]

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = []
        self.assertEquals(expected, actual_multalts)

    def test_build_coordinates_multAlts(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "T", "id=2")
        fileBrec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileArec1, fileArec2]),
                        MockVcfReader(records=[fileBrec1, fileBrec2])]

        actual_coordinates = merge2._build_coordinates(mock_readers)[0]

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = [fileArec2, fileBrec1]
        self.assertEquals(expected, actual_multalts)

    def test_build_merged_record_preserveSampleNamesAndOrder(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SD": {"foo":"bar1"},
                       "SC": {"foo":"bar2"}})
        samples2 = OD({"SB": {"foo":"bar3"},
                       "SA": {"foo":"bar4"}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        sample_list = ["SD", "SA", "SC", "SB"]
        actual_record = merge2._build_merged_record(coordinate, [record1, record2], sample_list)

        self.assertEquals(sample_list,  actual_record.sample_tag_values.keys())

    def test_build_merged_record_fillsMissingSamples(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SA": {"foo":"bar3"},
                       "SB": {"foo":"bar4"}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)

        sample_list = ["SA", "SB", "SC", "SD"]
        actual_record = merge2._build_merged_record(coordinate, [record1], sample_list)

        self.assertEquals(sample_list,  actual_record.sample_tag_values.keys())

    def test_build_merged_record_baseInfoCopiedFromCoordinate(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SA": {},
                       "SB": {}})
        samples2 = OD({"SC": {},
                       "SD": {}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        actual_record = merge2._build_merged_record(coordinate, [record1, record2], [])

        self.assertEquals("chr1", actual_record.chrom)
        self.assertEquals("1", actual_record.pos)
        self.assertEquals("A", actual_record.ref)
        self.assertEquals("C", actual_record.alt)
        self.assertEquals("baseInfo", actual_record.info)

    
    def test_build_merged_record_tagsOrdered(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        sampleA_tag_values = OD({"foo":"A1", "bar":"A2"})
        sampleB_tag_values = OD({"foo":"B1", "bar":"B2"})
        sampleC_tag_values = OD({"foo":"C1", "bar":"C2"})
        sampleD_tag_values = OD({"foo":"D1", "bar":"D2"})
        samples1 = OD({"SA": sampleA_tag_values,
                       "SB": sampleB_tag_values})
        samples2 = OD({"SC": sampleC_tag_values,
                       "SD": sampleD_tag_values})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        actual_record = merge2._build_merged_record(coordinate, [record1, record2], ["SA", "SB", "SC", "SD"])

        self.assertEquals(OD([("bar", "A2"), ("foo", "A1")]), actual_record.sample_tag_values["SA"])
        self.assertEquals(OD([("bar", "B2"), ("foo", "B1")]), actual_record.sample_tag_values["SB"])
        self.assertEquals(OD([("bar", "C2"), ("foo", "C1")]), actual_record.sample_tag_values["SC"])
        self.assertEquals(OD([("bar", "D2"), ("foo", "D1")]), actual_record.sample_tag_values["SD"])

    def test_build_merged_record_heterogeneousTags(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        sampleA_tag_values = OD({"foo":"A1", "bar":"A2"})
        sampleB_tag_values = OD({"foo":"B1", "bar":"B2"})
        sampleC_tag_values = OD({"baz":"C1"})
        sampleD_tag_values = OD({"baz":"D1"})
        samples1 = OD({"SA": sampleA_tag_values,
                       "SB": sampleB_tag_values})
        samples2 = OD({"SC": sampleC_tag_values,
                       "SD": sampleD_tag_values})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        actual_record = merge2._build_merged_record(coordinate, [record1, record2], ["SA", "SB", "SC", "SD"])

        self.assertEquals("chr1", actual_record.chrom)
        self.assertEquals("1", actual_record.pos)
        self.assertEquals("A", actual_record.ref)
        self.assertEquals("C", actual_record.alt)
        self.assertEquals("baseInfo", actual_record.info)
        self.assertEquals(set(["SA", "SB", "SC", "SD"]), set(actual_record.sample_tag_values.keys()))
        self.assertEquals(OD([("bar", "A2"), ("baz", "."), ("foo", "A1")]), actual_record.sample_tag_values["SA"])
        self.assertEquals(OD([("bar", "B2"), ("baz", "."), ("foo", "B1")]), actual_record.sample_tag_values["SB"])
        self.assertEquals(OD([("bar", "."), ("baz", "C1"), ("foo", ".")]), actual_record.sample_tag_values["SC"])
        self.assertEquals(OD([("bar", "."), ("baz", "D1"), ("foo", ".")]), actual_record.sample_tag_values["SD"])


#     def test_get_sample_tag_values(self):
#         OD = OrderedDict
#         sampleA_tag_values = OD({"foo":"A1", "bar":"A2"})
#         sampleB_tag_values = OD({"foo":"B1", "bar":"B2"})
#         samples1 = OD({"SA": sampleA_tag_values,
#                        "SB": sampleB_tag_values})
#         sampleC_tag_values = OD({"foo":"C1", "bar":"C2"})
#         sampleD_tag_values = OD({"foo":"D1", "bar":"D2"})
#         samples2 = OD({"SC": sampleC_tag_values,
#                        "SD": sampleD_tag_values})
#         record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
#         record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)
#         record3 = None
#         reader1 = MockBufferedReader([record1])
#         reader2 = MockBufferedReader([record2])
#         reader3 = MockBufferedReader([record3])
#         buffered_readers = [reader1, reader2, reader3]
#         merged_record = VcfRecord("chr1", "1", "A", "C")
# 
#         actual_sample_tag_values = merge2._get_tag_sample_values(buffered_readers, merged_record)
# 
#         self.assertEqual(set(["foo", "bar"]), set(actual_sample_tag_values.keys()))
#         self.assertEqual(OD({"SA":"A1","SB":"B1","SC":"C1","SD":"D1"}), actual_sample_tag_values["foo"])
#         self.assertEqual(OD({"SA":"A2","SB":"B2","SC":"C2","SD":"D2"}), actual_sample_tag_values["bar"])

    def Xtest_merge_records(self):
        coordinates = [VcfRecord("chrom","pos","ref","alt")]
        OD = OrderedDict 
        record1 = VcfRecord("chrom","pos","ref","alt", sample_tag_values=OD({"SA": OD({"foo":"A"}), "SB":OD({"foo":"B"})}))
        record2 = VcfRecord("chrom","pos","ref","alt", sample_tag_values=OD({"SC": OD({"foo":"C"}), "SD":OD({"foo":"D"})}))
        buffered_readers = [MockBufferedReader([record1]),MockBufferedReader([record2])]
        writer = MockFileWriter()
        merge2._merge_records(coordinates, buffered_readers, writer)
        self.assertEqual("chrom\tpos\tref\talt\t.\t.\t.\tfoo\tA\tB\tC\tD",writer.written[0])

    def Xtest_pull_matching_records(self):
        coordinate = VcfRecord("chrom","pos","ref","alt")
        OD = OrderedDict 
        record1 = VcfRecord("chrom","pos","ref","alt", sample_tag_values=OD({"SA": OD({"foo":"A"}), "SB":OD({"foo":"B"})}))
        record2 = VcfRecord("chrom","pos","ref","alt", sample_tag_values=OD({"SC": OD({"foo":"C"}), "SD":OD({"foo":"D"})}))
        record3 = VcfRecord("chrom","foo","ref","alt", sample_tag_values=OD({"SC": OD({"foo":"C"}), "SD":OD({"foo":"D"})}))
        buffered_readers = [MockBufferedReader([record1]),MockBufferedReader([record2]),MockBufferedReader([record3])]
        vcf_records = merge2._pull_matching_records(coordinate, buffered_readers)
        self.assertEqual([record1, record2], vcf_records)

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

class Merge2FunctionalTestCase(test_case.JacquardBaseTestCase):
    def Xtest_merge2(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "04_merge2")
            input_dir = os.path.join(module_testdir, "input")
            output_file = os.path.join(output_dir.path, "merged.vcf")

            command = ["merge2", input_dir, output_file, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)


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
        
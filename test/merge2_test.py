#pylint: disable=missing-docstring,line-too-long,too-many-public-methods
#pylint: disable=too-few-public-methods,too-many-instance-attributes
#pylint: disable=too-many-arguments,invalid-name,protected-access
from argparse import Namespace
from collections import OrderedDict
import os
from testfixtures import TempDirectory
import unittest

import jacquard.merge2 as merge2
import jacquard.vcf as vcf
import test_case as test_case
from jacquard.vcf import VcfRecord
import jacquard.utils as utils

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

        actual_coordinates = merge2._build_coordinates(mock_readers)

        expected = [fileArec1, fileArec2, fileBrec2]
        self.assertEquals(expected, actual_coordinates)

    def xtest_build_coordinates_sortsSampleNames(self):
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

        actual_coordinates = merge2._build_coordinates(mock_readers)

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

        actual_coordinates = merge2._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = [fileArec2, fileBrec1]
        self.assertEquals(expected, actual_multalts)

    def test_build_merged_record_onlyKeepJQTags(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SD": {"JQ_foo":"bar1",
                              "blahJQ_": "bar2"},
                       "SC": {"JQ_foo":"bar3",
                              "blah":"bar4"}})
        samples2 = OD({"SB": {"JQ_foo":"bar5"},
                       "SA": {"JQ_foo":"bar6"}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        sample_list = ["SA", "SB", "SC", "SD"]
        tags_to_keep = ["JQ_*"]
        actual_record = merge2._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

        self.assertEquals(OD([("JQ_foo", "bar6")]), actual_record.sample_tag_values["SA"])
        self.assertEquals(OD([("JQ_foo", "bar5")]), actual_record.sample_tag_values["SB"])
        self.assertEquals(OD([("JQ_foo", "bar3")]), actual_record.sample_tag_values["SC"])
        self.assertEquals(OD([("JQ_foo", "bar1")]), actual_record.sample_tag_values["SD"])

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
        tags_to_keep = ["foo"]
        actual_record = merge2._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

        self.assertEquals(sample_list,  actual_record.sample_tag_values.keys())

    def test_build_merged_record_fillsMissingSamples(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SA": {"foo":"bar3"},
                       "SB": {"foo":"bar4"}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)

        sample_list = ["SA", "SB", "SC", "SD"]
        tags_to_keep = ["foo"]
        actual_record = merge2._build_merged_record(coordinate, [record1], sample_list, tags_to_keep)

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

        actual_record = merge2._build_merged_record(coordinate, [record1, record2], [], [])

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

        sample_list = ["SA", "SB", "SC", "SD"]
        tags_to_keep = ["foo", "bar"]
        actual_record = merge2._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

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

        sample_list = ["SA", "SB", "SC", "SD"]
        tags_to_keep = ["foo", "bar", "baz"]
        actual_record = merge2._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

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

    def test_merge_records(self):
        coordinates = [VcfRecord("chrom","pos","ref","alt")]
        OD = OrderedDict
        record1 = VcfRecord("chrom","pos","ref","alt", sample_tag_values=OD({"SA": OD({"foo":"A"}), "SB":OD({"foo":"B"})}))
        record2 = VcfRecord("chrom","pos","ref","alt", sample_tag_values=OD({"SC": OD({"foo":"C"}), "SD":OD({"foo":"D"})}))
        buffered_readers = [MockBufferedReader([record1]),MockBufferedReader([record2])]
        writer = MockFileWriter()

        merge2._merge_records(coordinates, buffered_readers, writer, ["SA", "SB", "SC", "SD"], ["foo"])
        self.assertEqual("chrom\tpos\t.\tref\talt\t.\t.\t.\tfoo\tA\tB\tC\tD\n",writer.written[0])

    def test_pull_matching_records(self):
        coordinate = VcfRecord("chrom","pos","ref","alt")
        OD = OrderedDict
        record1 = VcfRecord("chrom","pos","ref","alt", sample_tag_values=OD({"SA": OD({"foo":"A"}), "SB":OD({"foo":"B"})}))
        record2 = VcfRecord("chrom","pos","ref","alt", sample_tag_values=OD({"SC": OD({"foo":"C"}), "SD":OD({"foo":"D"})}))
        buffered_readers = [MockBufferedReader([record1]), MockBufferedReader([record2])]

        vcf_records = merge2._pull_matching_records(coordinate, buffered_readers)
        self.assertEqual([record1, record2], vcf_records)

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

    def test_process_inputs(self):
        with TempDirectory() as input_dir:
            fileA = input_dir.write("fileA.vcf",
                                    "##source=strelka\n"
                                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_A\tSample_B\n")
            fileB = input_dir.write("fileB.vcf",
                                    "##source=strelka\n"
                                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_C\tSample_D\n")
            input_files = [fileA, fileB]
            actual_headers, actual_all_sample_names = merge2._process_inputs(input_files)

            expected_headers = ['##source=strelka',
                                "##jacquard.merge.file1=fileA.vcf(['Sample_A', 'Sample_B'])",
                                '##INFO=<ID=JQ_MULT_ALT_LOCUS,Number=0,Type=Flag,Description="dbSNP Membership",Source="Jacquard",Version="{}">'.format(utils.__version__),
                                "##jacquard.merge.file2=fileB.vcf(['Sample_C', 'Sample_D'])",
                                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfileA.vcf|Sample_A\tfileA.vcf|Sample_B\tfileB.vcf|Sample_C\tfileB.vcf|Sample_D"]

            expected_all_sample_names = ["fileA.vcf|Sample_A",
                                         "fileA.vcf|Sample_B",
                                         "fileB.vcf|Sample_C",
                                         "fileB.vcf|Sample_D"]

            self.assertEquals(5, len(expected_headers))
            self.assertEquals(expected_headers, actual_headers)
            self.assertEquals(4, len(actual_all_sample_names))
            self.assertEquals(expected_all_sample_names, actual_all_sample_names)

    def test_write_metaheaders(self):
        mock_writer = MockFileWriter()
        headers = ["##foo", "##bar", "#CHROM\tPOS"]
        exectution_context = ["##execution_context"]
        merge2._write_metaheaders(mock_writer, headers, exectution_context)

        self.assertEquals(["##foo\n##bar\n##execution_context\n", "#CHROM\tPOS\n"], mock_writer.written)

    def test_create_reader_lists(self):
        with TempDirectory() as input_dir:
            fileA = input_dir.write("fileA.vcf",
                                    "##source=strelka\n"
                                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_A\tSample_B\n"
                                    "chr1\t31\t.\tA\tT\t.\t.\t.\tDP\t23\t52\n")
            fileB = input_dir.write("fileB.vcf",
                                    "##source=strelka\n"
                                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_C\tSample_D\n"
                                    "chr2\t32\t.\tA\tT\t.\t.\t.\tDP\t24\t53\n")
            input_files = [fileA, fileB]
            buffered_readers, vcf_readers = merge2._create_reader_lists(input_files)

            for vcf_reader in vcf_readers:
                vcf_reader.close()

            self.assertEquals(2, len(buffered_readers))
            self.assertEquals(2, len(vcf_readers))

    def xtest_execute(self):
        vcf_content1 = ('''##source=strelka
##file1
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleA|SampleB
chr1|1|.|A|C|.|.|INFO|FORMAT|A_1|B_1
chr2|1|.|A|C|.|.|INFO|FORMAT|A_2|B_2
''').replace('|', "\t")
        vcf_content2 = ('''##source=strelka
##file2
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleC|SampleD
chr1|10|.|A|C|.|.|INFO|FORMAT|C_1|D_1
chr2|10|.|A|C|.|.|INFO|FORMAT|C_2|D_2
''').replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("fileA.vcf", vcf_content1)
            input_dir.write("fileB.vcf", vcf_content2)
            args = Namespace(input=input_dir.path,
                             output=os.path.join(output_dir.path, "fileC.vcf"))

            merge2.execute(args, ["##extra_header1", "##extra_header2"])

            output_dir.check("fileC.vcf")
            with open(os.path.join(output_dir.path, "fileC.vcf")) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()

        expected_output_headers = ["##source=strelka\n",
                                   "##file1\n",
                                   "##jacquard.merge.file1=fileA.vcf(['SampleA', 'SampleB'])\n",
                                   '##INFO=<ID=JQ_MULT_ALT_LOCUS,Number=0,Type=Flag,Description="dbSNP Membership",Source="Jacquard",Version="0.21">\n',
                                   "##file2\n",
                                   "##jacquard.merge.file2=fileB.vcf(['SampleC', 'SampleD'])\n",
                                   "##extra_header1\n",
                                   "##extra_header2\n",
                                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfileA.vcf|SampleA\tfileA.vcf|SampleB\tfileB.vcf|SampleC\tfileB.vcf|SampleD\n"]

        self.assertEquals(13, len(actual_output_lines))
        self.assertEquals(expected_output_headers, actual_output_lines[0:9])
        self.assertEquals("chr1\t1\t.\tA\tC\t.\t.\t.\tFORMAT\tA_1\tB_1\t.\t.\n", actual_output_lines[9])
        self.assertEquals("chr1\t10\t.\tA\tC\t.\t.\t.\tFORMAT\t.\t.\tC_1\tD_1\n", actual_output_lines[10])
        self.assertEquals("chr2\t1\t.\tA\tC\t.\t.\t.\tFORMAT\tA_2\tB_2\t.\t.\n", actual_output_lines[11])
        self.assertEquals("chr2\t10\t.\tA\tC\t.\t.\t.\tFORMAT\t.\t.\tC_2\tD_2\n", actual_output_lines[12])

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
        
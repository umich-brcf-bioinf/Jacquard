#pylint: disable=missing-docstring,line-too-long,too-many-public-methods
#pylint: disable=too-few-public-methods,too-many-instance-attributes
#pylint: disable=too-many-arguments,invalid-name,protected-access,global-statement
from __future__ import absolute_import
from argparse import Namespace
from collections import OrderedDict
from jacquard.vcf import VcfRecord
from test.vcf_test import MockVcfReader
from testfixtures import TempDirectory
import jacquard.logger
import jacquard.merge as merge
import jacquard.utils as utils
import jacquard.vcf as vcf
import os
import test.mock_logger
import test.test_case as test_case


class MockBufferedReader(object):
    def __init__(self, vcf_records):
        self.vcf_records_iter = iter(vcf_records)

    def next_if_equals(self, dummy):
        return self.vcf_records_iter.next()

class MockFileWriter(object):
    def __init__(self):
        self.written = []

    def write(self, text):
        self.written.append(text)

class MergeTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(MergeTestCase, self).setUp()
        merge.logger = test.mock_logger

    def tearDown(self):
        test.mock_logger.reset()
        merge.logger = jacquard.logger
        super(MergeTestCase, self).tearDown()

    def test_validate_arguments(self):
        with TempDirectory() as input_file, TempDirectory() as output_file:
            input_file.write("fileA.vcf",
                             "##source=strelka\n"
                             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_A\tSample_B\n"
                             "chr1\t31\t.\tA\tT\t.\t.\t.\tDP\t23\t52\n")
            input_file.write("fileB.vcf",
                             "##source=strelka\n"
                             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_C\tSample_D\n"
                             "chr2\t32\t.\tA\tT\t.\t.\t.\tDP\t24\t53\n")

            args = Namespace(input=input_file.path,
                             output=os.path.join(output_file.path, "output.vcf"),
                             tags=None)
            actual_readers_to_writers, out_files, buffered_readers, format_tag_regex = merge._validate_arguments(args)

            for values in actual_readers_to_writers.values():
                for vcf_reader in values:
                    vcf_reader.close()

            self.assertEquals(1, len(actual_readers_to_writers))
            self.assertEquals(0, len(out_files))
            self.assertEquals(2, len(buffered_readers))
            self.assertEquals(merge._DEFAULT_INCLUDED_FORMAT_TAGS, format_tag_regex)
            self.assertRegexpMatches(actual_readers_to_writers.keys()[0].output_filepath, ".merged.vcf")

    def test_predict_output(self):
        with TempDirectory() as input_file, TempDirectory() as output_file:
            input_file.write("A.normalized.jacquardTags.HCsomatic.vcf", "##source=strelka\n#colHeader")
            input_file.write("B.normalized.jacquardTags.HCsomatic.vcf", "##source=strelka\n#colHeader")
            output_file.write("merged.vcf", "##source=strelka\n#colHeader")
            args = Namespace(input=input_file.path,
                             output=os.path.join(output_file.path, "merged.vcf"))

            desired_output_files = merge._predict_output(args)
            expected_desired_output_files = set(["merged.vcf"])

            self.assertEquals(expected_desired_output_files, desired_output_files)

    def test_build_coordinates(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "G", "id=2")
        fileBrec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileArec1, fileArec2]),
                        MockVcfReader(records=[fileBrec1, fileBrec2])]

        actual_coordinates = merge._build_coordinates(mock_readers)

        expected = [fileArec1, fileArec2, fileBrec2]
        self.assertEquals(expected, actual_coordinates)

    def test_build_coordinates_unsorted(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "5", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "16", "A", "G", "id=2")
        fileBrec2 = vcf.VcfRecord("chr2", "12", "G", "C")

        mock_readers = [MockVcfReader(input_filepath="fileA.vcf",
                                      records=[fileArec1, fileArec2]),
                        MockVcfReader(input_filepath="fileB.vcf",
                                      records=[fileBrec1, fileBrec2])]

        self.assertRaisesRegexp(utils.JQException,
                                "One or more VCF files were not sorted.*",
                                merge._build_coordinates, mock_readers)
        actual_log_errors = test.mock_logger.messages["ERROR"]
        self.assertRegexpMatches(actual_log_errors[0],
                                 r"VCF file:chrom:pos \[fileB.vcf:chr2:12\] is out of order")

    def test_build_coordinates_multAltsEmpty(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "G", "id=2")
        fileBrec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileArec1, fileArec2]),
                        MockVcfReader(records=[fileBrec1, fileBrec2])]

        actual_coordinates = merge._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        self.assertEquals([], actual_multalts)

    def test_build_coordinates_flagsMultAltsFromDistinctFiles(self):
        fileA_rec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileA_rec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileB_rec1 = vcf.VcfRecord("chr2", "12", "A", "T", "id=2")
        fileB_rec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileA_rec1, fileA_rec2]),
                        MockVcfReader(records=[fileB_rec1, fileB_rec2])]

        actual_coordinates = merge._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = [fileA_rec2, fileB_rec1]
        self.assertEquals(expected, actual_multalts)

    def test_build_coordinates_flagsMultAltsWithinFile(self):
        fileA_rec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileA_rec2 = vcf.VcfRecord("chr2", "12", "A", "G,T", "id=1")
        fileB_rec1 = vcf.VcfRecord("chr3", "12", "A", "T", "id=2")
        fileB_rec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileA_rec1, fileA_rec2]),
                        MockVcfReader(records=[fileB_rec1, fileB_rec2])]

        actual_coordinates = merge._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = [fileA_rec2]
        self.assertEquals(expected, actual_multalts)

    def test_build_coordinates_flagsMultAltsWithDistinctRefs(self):
        fileA_rec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileA_rec2 = vcf.VcfRecord("chr2", "2", "A", "G", "id=1")
        fileB_rec1 = vcf.VcfRecord("chr2", "2", "AT", "T", "id=2")
        fileB_rec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileA_rec1, fileA_rec2]),
                        MockVcfReader(records=[fileB_rec1, fileB_rec2])]

        actual_coordinates = merge._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = [fileA_rec2, fileB_rec1]
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
        tags_to_keep = ["JQ_foo", "JQ_foo"]
        actual_record = merge._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

        self.assertEquals(OD([("JQ_foo", "bar6")]), actual_record.sample_tag_values["SA"])
        self.assertEquals(OD([("JQ_foo", "bar5")]), actual_record.sample_tag_values["SB"])
        self.assertEquals(OD([("JQ_foo", "bar3")]), actual_record.sample_tag_values["SC"])
        self.assertEquals(OD([("JQ_foo", "bar1")]), actual_record.sample_tag_values["SD"])

    def test_build_merged_record_redundantPatientNames(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"PA|NORMAL": {"JQ_vs":"1"},
                       "PA|TUMOR": {"JQ_vs":"2"}})
        samples2 = OD({"PA|NORMAL": {"JQ_mt":"3"},
                       "PA|TUMOR": {"JQ_mt":"4"}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        sample_list = ["PA|NORMAL", "PA|TUMOR"]
        tags_to_keep = ["JQ_vs", "JQ_mt"]
        actual_record = merge._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

        self.assertEquals(OD([("JQ_mt", "3"), ("JQ_vs", "1")]), actual_record.sample_tag_values["PA|NORMAL"])
        self.assertEquals(OD([("JQ_mt", "4"), ("JQ_vs", "2")]), actual_record.sample_tag_values["PA|TUMOR"])

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
        actual_record = merge._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

        self.assertEquals(sample_list, actual_record.sample_tag_values.keys())

    def test_build_merged_record_fillsMissingSamples(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SA": {"foo":"bar3"},
                       "SB": {"foo":"bar4"}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)

        sample_list = ["SA", "SB", "SC", "SD"]
        tags_to_keep = ["foo"]
        actual_record = merge._build_merged_record(coordinate, [record1], sample_list, tags_to_keep)

        self.assertEquals(sample_list, actual_record.sample_tag_values.keys())

    def test_build_merged_record_baseInfoCopiedFromCoordinate(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SA": {},
                       "SB": {}})
        samples2 = OD({"SC": {},
                       "SD": {}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        actual_record = merge._build_merged_record(coordinate, [record1, record2], [], [])

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
        actual_record = merge._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

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
        actual_record = merge._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

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

    def test_merge_records(self):
        coordinates = [VcfRecord("chrom", "pos", "ref", "alt")]
        OD = OrderedDict
        record1 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=OD({"SA": OD({"foo":"A"}), "SB":OD({"foo":"B"})}))
        record2 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=OD({"SC": OD({"foo":"C"}), "SD":OD({"foo":"D"})}))
        buffered_readers = [MockBufferedReader([record1]), MockBufferedReader([record2])]
        writer = MockFileWriter()

        merge._merge_records(coordinates, buffered_readers, writer, ["SA", "SB", "SC", "SD"], ["foo"])
        self.assertEqual("chrom\tpos\t.\tref\talt\t.\t.\t.\tfoo\tA\tB\tC\tD\n", writer.written[0])

    def test_pull_matching_records(self):
        coordinate = VcfRecord("chrom", "pos", "ref", "alt")
        OD = OrderedDict
        record1 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=OD({"SA": OD({"foo":"A"}), "SB":OD({"foo":"B"})}))
        record2 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=OD({"SC": OD({"foo":"C"}), "SD":OD({"foo":"D"})}))
        buffered_readers = [MockBufferedReader([record1]), MockBufferedReader([record2])]

        vcf_records = merge._pull_matching_records(coordinate, buffered_readers)
        self.assertEqual([record1, record2], vcf_records)

    def test_build_sample_list_simpleSampleList(self):
        reader1 = MockVcfReader("PA.foo.vcf",
                                sample_names=["Sample_A", "Sample_B"])
        reader2 = MockVcfReader("PA.bar.vcf",
                                sample_names=["Sample_A", "Sample_B"])
        reader3 = MockVcfReader("PB.vcf",
                                sample_names=["Sample_C", "Sample_D"])
        readers = [reader1, reader2, reader3]
        actual_sample_names, dummy = merge._build_sample_list(readers)

        expected_sample_names = ["PA|Sample_A",
                                 "PA|Sample_B",
                                 "PB|Sample_C",
                                 "PB|Sample_D"]
        self.assertEquals(expected_sample_names, actual_sample_names)

    def test_build_sample_list_patientNamesNaturalOrdered(self):
        reader1 = MockVcfReader("P10.foo.vcf", sample_names=["S10", "S2"])
        reader2 = MockVcfReader("P1A.bar.vcf", sample_names=["S10", "S2"])
        reader3 = MockVcfReader("P1.bar.vcf", sample_names=["S10", "S2"])
        reader4 = MockVcfReader("P2.vcf", sample_names=["S10", "S2"])
        readers = [reader1, reader2, reader3, reader4]
        actual_sample_names, dummy = merge._build_sample_list(readers)

        expected_sample_names = ["P1|S2", "P1|S10",
                                 "P1A|S2", "P1A|S10",
                                 "P2|S2", "P2|S10",
                                 "P10|S2", "P10|S10"]
        self.assertEquals(expected_sample_names, actual_sample_names)

    def test_build_sample_list_sampleNamesOrdered(self):
        reader1 = MockVcfReader("P10.foo.vcf", sample_names=["S10", "S1", "S100", "S2"])
        reader2 = MockVcfReader("P1.foo.vcf", sample_names=["S0", "S1"])
        readers = [reader1, reader2]
        actual_sample_names, dummy = merge._build_sample_list(readers)

        expected_sample_names = ["P1|S0", "P1|S1",
                                 "P10|S1", "P10|S2", "P10|S10", "P10|S100"]
        self.assertEquals(expected_sample_names, actual_sample_names)

    def test_build_sample_list_mergeMetaheaders(self):
        reader1 = MockVcfReader("PA.foo.vcf",
                                sample_names=["Sample_A", "Sample_B"])
        reader2 = MockVcfReader("PA.bar.vcf",
                                sample_names=["Sample_A", "Sample_B"])
        reader3 = MockVcfReader("PB.vcf",
                                sample_names=["Sample_C", "Sample_D"])
        readers = [reader1, reader2, reader3]

        dummy, actual_metaheaders = merge._build_sample_list(readers)

        expected_metaheaders = ["##jacquard.merge.sample=<Column=1,Name=PA|Sample_A,Source=PA.foo.vcf|PA.bar.vcf>",
                                "##jacquard.merge.sample=<Column=2,Name=PA|Sample_B,Source=PA.foo.vcf|PA.bar.vcf>",
                                "##jacquard.merge.sample=<Column=3,Name=PB|Sample_C,Source=PB.vcf>",
                                "##jacquard.merge.sample=<Column=4,Name=PB|Sample_D,Source=PB.vcf>"]
        self.assertEquals(4, len(actual_metaheaders))
        self.assertEquals(expected_metaheaders, actual_metaheaders)

    def test_create_reader_lists(self):
        with TempDirectory() as input_file:
            fileA = input_file.write("fileA.vcf",
                                     "##source=strelka\n"
                                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_A\tSample_B\n"
                                     "chr1\t31\t.\tA\tT\t.\t.\t.\tDP\t23\t52\n")
            fileB = input_file.write("fileB.vcf",
                                     "##source=strelka\n"
                                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_C\tSample_D\n"
                                     "chr2\t32\t.\tA\tT\t.\t.\t.\tDP\t24\t53\n")
            input_files = [fileA, fileB]
            buffered_readers, vcf_readers = merge._create_reader_lists(input_files)

            for vcf_reader in vcf_readers:
                vcf_reader.close()

            self.assertEquals(2, len(buffered_readers))
            self.assertEquals(2, len(vcf_readers))

    def test_build_info_tags_sorts(self):
        records = [VcfRecord("1", "42", "A", "C", info="foo"),
                   VcfRecord("1", "43", "A", "C", info="bar")]
        actual_tags = merge._build_info_tags(records)
        self.assertEquals(["bar", "foo"], actual_tags)

    def test_build_info_tags_noDuplicates(self):
        records = [VcfRecord("1", "42", "A", "C", info="foo"),
                   VcfRecord("1", "43", "A", "C", info="bar;foo")]
        actual_tags = merge._build_info_tags(records)
        self.assertEquals(["bar", "foo"], actual_tags)

    def test_build_info_tags_empty(self):
        records = [VcfRecord("1", "42", "A", "C", info="")]
        actual_tags = merge._build_info_tags(records)
        self.assertEquals([], actual_tags)

    def test_build_info_tags_null(self):
        records = [VcfRecord("1", "42", "A", "C", info=".")]
        actual_tags = merge._build_info_tags(records)
        self.assertEquals([], actual_tags)

    def test_build_format_tags_to_keep_filtersAndSorts(self):
        meta_headers = ['##FORMAT=<ID=JQ1>',
                        '##FORMAT=<ID=JQ2>',
                        '##FORMAT=<ID=JQ3>']
        vcf_reader = MockVcfReader(metaheaders=meta_headers)
        format_tags = merge._build_format_tags(["JQ[12]"], [vcf_reader])

        self.assertEquals(["JQ1", "JQ2"], format_tags)

    def test_build_format_tags_to_keep_uniqueSetAcrossMultipleReaders(self):
        reader1 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1>'])
        reader2 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1>', '##FORMAT=<ID=JQ2>'])
        reader3 = MockVcfReader(metaheaders=['##INFO=<ID=JQ3>'])

        format_tags = merge._build_format_tags(["JQ[123]"], [reader1, reader2, reader3])

        self.assertEquals(["JQ1", "JQ2"], format_tags)


    def test_build_contigs(self):
        records = [VcfRecord("1", "42", "A", "C"),
                   VcfRecord("4", "42", "A", "C"),
                   VcfRecord("4", "44", "A", "C"),
                   VcfRecord("15", "42", "A", "C")]
        actual_contigs = merge._build_contigs(records)

        expected_contigs = ["1", "4", "15"]
        self.assertEquals(expected_contigs, actual_contigs)

    def test_compile_metaheaders_preservesExistingMetaheaders(self):
        reader1 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1>'])
        vcf_readers = [reader1]
        all_sample_names = []
        contigs_to_keep = []
        format_tags_to_keep = []
        info_tags_to_keep = []
        existing_metaheaders = ["##existing1", "##existing2"]
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)
        self.assertEquals(3, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertEquals("##existing1", metaheaders_iter.next())
        self.assertEquals("##existing2", metaheaders_iter.next())
        self.assertRegexpMatches(metaheaders_iter.next(), "^#CHROM")

    def test_compile_metaheaders_retainsFormatMetaheaders(self):
        reader1 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1,Description="JQA">',
                                             '##FORMAT=<ID=JQ2,Description="JQB">',
                                             '##FORMAT=<ID=JQ3,Description="JQC">'])
        vcf_readers = [reader1]
        all_sample_names = []
        contigs_to_keep = []
        format_tags_to_keep = ["JQ1", "JQ3"]
        info_tags_to_keep = []
        existing_metaheaders = []
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)
        self.assertEquals(3, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(metaheaders_iter.next(), "##FORMAT.*JQA")
        self.assertRegexpMatches(metaheaders_iter.next(), "##FORMAT.*JQC")
        self.assertRegexpMatches(metaheaders_iter.next(), "^#CHROM")

    def test_compile_metaheaders_retainsInfoMetaheaders(self):
        reader1 = MockVcfReader(metaheaders=['##INFO=<ID=JQ1,Description="JQA">',
                                             '##INFO=<ID=JQ2,Description="JQB">',
                                             '##INFO=<ID=JQ3,Description="JQC">'])
        vcf_readers = [reader1]
        all_sample_names = []
        contigs_to_keep = []
        format_tags_to_keep = []
        info_tags_to_keep = ["JQ1", "JQ3", merge._MULT_ALT_TAG]
        existing_metaheaders = []
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)
        self.assertEquals(4, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(metaheaders_iter.next(), "##INFO.*JQA")
        self.assertRegexpMatches(metaheaders_iter.next(), "##INFO.*JQC")
        self.assertRegexpMatches(metaheaders_iter.next(), merge._MULT_ALT_HEADER)
        self.assertRegexpMatches(metaheaders_iter.next(), "^#CHROM")

    def test_compile_metaheaders_retainsContigMetaheaders(self):
        reader1 = MockVcfReader(metaheaders=['##contig=<ID=chr1,Description="chromosome 1">',
                                             '##contig=<ID=chr2,Description="chromosome 2">',
                                             '##contig=<ID=chr4,Description="chromosome 2">',
                                             '##contig=<ID=chr11,Description="chromosome 11">'])
        vcf_readers = [reader1]
        all_sample_names = []
        contigs_to_keep = ["chr1", "chr2", "chr11"]
        format_tags_to_keep = []
        info_tags_to_keep = []
        existing_metaheaders = []
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)

        self.assertEquals(4, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(metaheaders_iter.next(), "##contig.*chromosome 1")
        self.assertRegexpMatches(metaheaders_iter.next(), "##contig.*chromosome 2")
        self.assertRegexpMatches(metaheaders_iter.next(), "##contig.*chromosome 11")
        self.assertRegexpMatches(metaheaders_iter.next(), "^#CHROM")

    def test_compile_metaheaders_noContigMetaheaders(self):
        reader1 = MockVcfReader(metaheaders=['##foo'])
        vcf_readers = [reader1]
        all_sample_names = []
        contigs_to_keep = ["chr1", "chr2", "chr11"]
        format_tags_to_keep = []
        info_tags_to_keep = []
        existing_metaheaders = []
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)

        self.assertEquals(1, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(metaheaders_iter.next(), "^#CHROM")

    def test_compile_metaheaders_addsSampleNamesToColumnHeader(self):
        reader1 = MockVcfReader()
        vcf_readers = [reader1]
        all_sample_names = ["SA", "SB"]
        contigs_to_keep = []
        format_tags_to_keep = []
        info_tags_to_keep = []
        existing_metaheaders = []
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)
        self.assertEquals(1, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(metaheaders_iter.next(), "#CHROM.*SA\tSB")

    def test_compile_metaheaders_ordersHeaders(self):
        reader1 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1,Description="JQA">'])
        reader2 = MockVcfReader(metaheaders=['##INFO=<ID=JQ2,Description="JQB">'])

        vcf_readers = [reader1, reader2]
        all_sample_names = []
        contigs_to_keep = []
        format_tags_to_keep = ["JQ1"]
        info_tags_to_keep = ["JQ2"]
        incoming_metaheaders = ["##foo"]
        actual_metaheaders = merge._compile_metaheaders(incoming_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)
        self.assertEquals(4, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(metaheaders_iter.next(), "##foo")
        self.assertRegexpMatches(metaheaders_iter.next(), "##INFO.*JQB")
        self.assertRegexpMatches(metaheaders_iter.next(), "##FORMAT.*JQA")
        self.assertRegexpMatches(metaheaders_iter.next(), "#CHROM.*")

    def test_execute(self):
        vcf_content1 = ('''##source=strelka
##contig=<ID=chr1,Number=1>
##contig=<ID=chr2,Number=1>
##FORMAT=<ID=JQ_Foo1,Number=1,Type=Float,Description="foo",Source="Jacquard",Version=0.21>
##FORMAT=<ID=JQ_Bar1,Number=1,Type=Float,Description="bar",Source="Jacquard",Version=0.21>
##file1
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleA|SampleB
chr1|1|.|A|C|.|.|INFO|JQ_Foo1:JQ_Bar1|A_1_1:A_1_2|B_1_1:B_1_2
chr1|1|.|A|T|.|.|INFO|JQ_Foo1|A_2|B_2
chr2|1|.|A|C|.|.|INFO|JQ_Foo1:JQ_Bar1|A_3_1:A_3_2|B_3_1:B_3_2
''').replace('|', "\t")
        vcf_content2 = ('''##source=strelka
##contig=<ID=chr1,Number=1>
##contig=<ID=chr2,Number=1>
##file2
##FORMAT=<ID=JQ_Foo2,Number=1,Type=Float,Description="foo",Source="Jacquard",Version=0.21>
##FORMAT=<ID=JQ_Bar2,Number=1,Type=Float,Description="bar",Source="Jacquard",Version=0.21>
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleA|SampleB
chr1|10|.|A|C|.|.|INFO|JQ_Foo2|C_1_1|D_1_2
chr2|10|.|A|C|.|.|INFO|JQ_Bar2|C_2|D_2
''').replace('|', "\t")

        with TempDirectory() as input_file, TempDirectory() as output_file:
            input_file.write("P1.fileA.vcf", vcf_content1)
            input_file.write("P1.fileB.vcf", vcf_content2)
            args = Namespace(input=input_file.path,
                             output=os.path.join(output_file.path, "fileC.vcf"),
                             tags=None)

            merge.execute(args, ["##execution_header1", "##execution_header2"])

            output_file.check("fileC.vcf")
            with open(os.path.join(output_file.path, "fileC.vcf")) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()

        self.assertEquals(18, len(actual_output_lines))
        actual_lines_iter = iter(actual_output_lines)
        self.assertEquals("##fileformat=VCFv4.2\n", actual_lines_iter.next())
        self.assertEquals("##execution_header1\n", actual_lines_iter.next())
        self.assertEquals("##execution_header2\n", actual_lines_iter.next())
        self.assertRegexpMatches(actual_lines_iter.next(), "##jacquard.merge.sample=<Column=1.*>\n")
        self.assertRegexpMatches(actual_lines_iter.next(), "##jacquard.merge.sample=<Column=2.*>\n")
        self.assertRegexpMatches(actual_lines_iter.next(), "##contig=<ID=chr1,Number=1>\n")
        self.assertRegexpMatches(actual_lines_iter.next(), "##contig=<ID=chr2,Number=1>\n")
        self.assertRegexpMatches(actual_lines_iter.next(), '##INFO=<ID=JQ_MULT_ALT_LOCUS.*>\n')
        self.assertRegexpMatches(actual_lines_iter.next(), '##FORMAT=<ID=JQ_Bar1.*>\n')
        self.assertRegexpMatches(actual_lines_iter.next(), '##FORMAT=<ID=JQ_Bar2.*>\n')
        self.assertRegexpMatches(actual_lines_iter.next(), '##FORMAT=<ID=JQ_Foo1.*>\n')
        self.assertRegexpMatches(actual_lines_iter.next(), '##FORMAT=<ID=JQ_Foo2.*>\n')
        self.assertRegexpMatches(actual_lines_iter.next(), "#CHROM\t.*P1|SampleA\tP1|SampleB\n")
        self.assertEquals("chr1\t1\t.\tA\tC\t.\t.\tJQ_MULT_ALT_LOCUS\tJQ_Bar1:JQ_Foo1\tA_1_2:A_1_1\tB_1_2:B_1_1\n", actual_lines_iter.next())
        self.assertEquals("chr1\t1\t.\tA\tT\t.\t.\tJQ_MULT_ALT_LOCUS\tJQ_Foo1\tA_2\tB_2\n", actual_lines_iter.next())
        self.assertEquals("chr1\t10\t.\tA\tC\t.\t.\t.\tJQ_Foo2\tC_1_1\tD_1_2\n", actual_lines_iter.next())
        self.assertEquals("chr2\t1\t.\tA\tC\t.\t.\t.\tJQ_Bar1:JQ_Foo1\tA_3_2:A_3_1\tB_3_2:B_3_1\n", actual_lines_iter.next())
        self.assertEquals("chr2\t10\t.\tA\tC\t.\t.\t.\tJQ_Bar2\tC_2\tD_2\n", actual_lines_iter.next())

    def test_execute_includeFormatIds(self):
        vcf_content1 = ('''##source=strelka
##contig=<ID=chr1,Number=1>
##contig=<ID=chr2,Number=1>
##FORMAT=<ID=JQ_Foo,Number=1,Type=Float,Description="foo",Source="Jacquard",Version=0.21>
##FORMAT=<ID=JQ_Foo1,Number=1,Type=Float,Description="foo",Source="Jacquard",Version=0.21>
##FORMAT=<ID=Bar,Number=1,Type=Float,Description="bar",Source="Jacquard",Version=0.21>
##file1
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleA|SampleB
chr1|1|.|A|C|.|.|INFO|JQ_Foo1:JQ_Bar1|A_1_1:A_1_2|B_1_1:B_1_2
chr1|1|.|A|T|.|.|INFO|JQ_Foo1|A_2|B_2
chr2|1|.|A|C|.|.|INFO|JQ_Foo1:JQ_Bar1|A_3_1:A_3_2|B_3_1:B_3_2
''').replace('|', "\t")
        with TempDirectory() as input_file, TempDirectory() as output_file:
            input_file.write("P1.fileA.vcf", vcf_content1)
            args = Namespace(input=input_file.path,
                             output=os.path.join(output_file.path, "fileB.vcf"),
                             tags="JQ_Foo,Bar")

            merge.execute(args, ["##extra_header1", "##extra_header2"])

            output_file.check("fileB.vcf")
            with open(os.path.join(output_file.path, "fileB.vcf")) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()

        expected_output_headers = ["##fileformat=VCFv4.2\n",
                                   "##extra_header1\n",
                                   "##extra_header2\n",
                                   "##jacquard.merge.sample=<Column=1,Name=P1|SampleA,Source=P1.fileA.vcf>\n",
                                   "##jacquard.merge.sample=<Column=2,Name=P1|SampleB,Source=P1.fileA.vcf>\n",
                                   '##contig=<ID=chr1,Number=1>\n',
                                   '##contig=<ID=chr2,Number=1>\n',
                                   '##INFO=<ID=JQ_MULT_ALT_LOCUS,Number=0,Type=Flag,Description="dbSNP Membership",Source="Jacquard",Version="0.21">\n',
                                   '##FORMAT=<ID=Bar,Number=1,Type=Float,Description="bar",Source="Jacquard",Version=0.21>\n',
                                   '##FORMAT=<ID=JQ_Foo,Number=1,Type=Float,Description="foo",Source="Jacquard",Version=0.21>\n',
                                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1|SampleA\tP1|SampleB\n"]

        self.assertEquals(expected_output_headers, actual_output_lines[0:len(expected_output_headers)])

    def test_sort_vcfs_orderedVcfsPassThrough(self):
        record1 = vcf.VcfRecord("chr1", "42", "A", "C")
        record2 = vcf.VcfRecord("chr2", "42", "A", "C")
        record3 = vcf.VcfRecord("chr3", "42", "A", "C")
        vcf_readerA = MockVcfReader(records=[record1, record2, record3])
        vcf_readerB = MockVcfReader(records=[record1, record2, record3])

        input_readers = [vcf_readerA, vcf_readerB]
        with TempDirectory() as temp_dir:
            actual_readers = merge._sort_readers(list(input_readers),
                                                 temp_dir.path)

        self.assertEquals(actual_readers, input_readers)

    def test_sort_readers_vcfsResortedAsNecessary(self):
        #pylint: disable=too-many-locals
        record1 = vcf.VcfRecord("chr1", "42", "A", "C")
        record2 = vcf.VcfRecord("chr2", "42", "A", "C")
        record3 = vcf.VcfRecord("chr3", "42", "A", "C")
        vcf_readerA = MockVcfReader(records=[record1, record2, record3])
        input_metaheaders = ["##foo", "##bar"]
        input_column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SX|SY")
        vcf_readerB = MockVcfReader(input_filepath="unsorted.vcf",
                                    metaheaders=input_metaheaders,
                                    column_header=input_column_header,
                                    records=[record3, record1, record2])

        input_readers = [vcf_readerA, vcf_readerB]
        with TempDirectory() as temp_dir:
            actual_readers = merge._sort_readers(list(input_readers),
                                                 temp_dir.path)

            self.assertNotEquals(actual_readers, input_readers)
            self.assertIn(vcf_readerA, actual_readers)
            self.assertNotIn(vcf_readerB, actual_readers)
            self.assertEquals(2, len(actual_readers))
            actual_readerB = actual_readers[1]
            actual_readerB.open()
            actual_metaheaders = list(actual_readerB.metaheaders)
            actual_column_header = actual_readerB.column_header
            actual_records = [vcf_record for vcf_record in actual_readerB.vcf_records()]
            actual_readerB.close()


        self.assertEquals(input_metaheaders, actual_metaheaders)
        self.assertEquals(input_column_header, actual_column_header)
        self.assertEquals(record1.asText(), actual_records[0].asText())
        self.assertEquals(record2.asText(), actual_records[1].asText())
        self.assertEquals(record3.asText(), actual_records[2].asText())
        actual_log_infos = test.mock_logger.messages["INFO"]
        self.assertEquals(1, len(actual_log_infos))
        self.assertRegexpMatches(actual_log_infos[0], r"Sorting vcf \[unsorted.vcf\]")

class MergeFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_merge(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "04_merge")
            input_file = os.path.join(module_testdir, "input")
            output_file = os.path.join(output_dir.path, "merged.vcf")

            command = ["merge", input_file, output_file, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

class BufferedReaderTestCase(test_case.JacquardBaseTestCase):
    def test_get_sample_info_advancesCurrentElementWhenMatched(self):
        reader = [1, 5, 10, 15]
        buffered_reader = merge._BufferedReader(iter(reader))

        self.assertEquals(None, buffered_reader.next_if_equals(0))
        self.assertEquals(None, buffered_reader.next_if_equals(5))
        self.assertEquals(1, buffered_reader.next_if_equals(1))
        self.assertEquals(None, buffered_reader.next_if_equals(1))
        self.assertEquals(5, buffered_reader.next_if_equals(5))
        self.assertEquals(10, buffered_reader.next_if_equals(10))
        self.assertEquals(15, buffered_reader.next_if_equals(15))
        self.assertEquals(None, buffered_reader.next_if_equals(15))
        self.assertEquals(None, buffered_reader.next_if_equals(42))
        
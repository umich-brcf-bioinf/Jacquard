# pylint: disable=line-too-long,too-many-public-methods,too-few-public-methods
# pylint: disable=invalid-name,global-statement,too-many-format-args
from jacquard import __version__
from jacquard.utils import JQException
from jacquard.vcf import VcfRecord
import jacquard.variant_callers.common_tags as common_tags
import jacquard.variant_callers.strelka as strelka
import os
import unittest

#TODO: (cgates): Lots of PEP8 cleanup in this class
ORIGINAL_REPORTED_TAG = None
ORIGINAL_PASSED_TAG = None

class MockCommonTag(object):
    def __init__(self, input_caller_name):
        self.input_caller_name = input_caller_name

class MockWriter(object):
    def __init__(self):
        self._content = []
        self.opened = False
        self.closed = False

    def open(self):
        self.opened = True

    def write(self, content):
        self._content.extend(content.splitlines())

    def lines(self):
        return self._content

    def close(self):
        self.closed = True

class MockFileReader(object):
    def __init__(self, input_filepath="/foo/mockFileReader.txt", content=None):
        if not content:
            content = []
        self.input_filepath = input_filepath
        self.file_name = os.path.basename(input_filepath)
        self._content = content
        self.open_was_called = False
        self.close_was_called = False

    def open(self):
        self.open_was_called = True

    def read_lines(self):
        for line in self._content:
            yield line

    def close(self):
        self.close_was_called = True

class CommonTagTestCase(unittest.TestCase):
    def setUp(self):
        global ORIGINAL_REPORTED_TAG
        global ORIGINAL_PASSED_TAG
        ORIGINAL_REPORTED_TAG = common_tags.ReportedTag
        ORIGINAL_PASSED_TAG = common_tags.PassedTag
        common_tags.ReportedTag = MockCommonTag
        common_tags.PassedTag = MockCommonTag

    def tearDown(self):
        common_tags.ReportedTag = ORIGINAL_REPORTED_TAG
        common_tags.PassedTag = ORIGINAL_PASSED_TAG

    def test_reported_tag(self):
        strelka_instance = strelka.Strelka()
        reported_tag = strelka_instance.tags[0]
        passed_tag = strelka_instance.tags[1]
        self.assertEquals("JQ_SK_", reported_tag.input_caller_name)
        self.assertEquals("JQ_SK_", passed_tag.input_caller_name)

class AlleleFreqTagTestCase(unittest.TestCase):

    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_SK_AF,Number=A,Type=Float,Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth. Uses (TIR tier 2)/DP2 if available, otherwise uses (ACGT tier2 depth) / DP2)",Source="Jacquard",Version=0.21>'.format(strelka.JQ_STRELKA_TAG, __version__), strelka._AlleleFreqTag().metaheader)

    def test_format_missingAFTag(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        originalVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())

    def test_format_AUTag(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|A,C|QUAL|FILTER|INFO|AU:CU:GU:TU|1,2:3,4:5,6:7,8|9,10:11,12:13,14:15,16\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|A,C|QUAL|FILTER|INFO|AU:CU:GU:TU:{0}AF|1,2:3,4:5,6:7,8:0.1,0.2|9,10:11,12:13,14:15,16:0.19,0.23\n".format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_AFTag_noAlt(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|.|QUAL|FILTER|INFO|AU:CU:GU:TU|1,2:3,4:5,6:7,8|9,10:11,12:13,14:15,16\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|.|QUAL|FILTER|INFO|AU:CU:GU:TU:{0}AF|1,2:3,4:5,6:7,8:.|9,10:11,12:13,14:15,16:.\n".format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_TIRTag(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:TIR|10:3,4|20:11,7\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:TIR:{0}AF|10:3,4:0.4|20:11,7:0.35\n".format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class DepthTagTestCase(unittest.TestCase):

    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}DP,Number=1,Type=Float,Description="Jacquard depth for Strelka (uses DP2 if available, otherwise uses ACGT tier2 depth)",Source="Jacquard",Version={1}>'.format(strelka.JQ_STRELKA_TAG, __version__), strelka._DepthTag().metaheader)

    def test_format_missingTag(self):
        tag = strelka._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        originalVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())

    def test_format_DP2Tag(self):
        tag = strelka._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:F2:F3|2:SA.2:SA.3|4:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:F2:F3:{0}DP|2:SA.2:SA.3:2|4:SB.2:SB.3:4\n".format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_AUTag(self):
        tag = strelka._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|AU:CU:GU:TU|1,2:3,4:5,6:7,8|9,10:11,12:13,14:15,16\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|AU:CU:GU:TU:{0}DP|1,2:3,4:5,6:7,8:20|9,10:11,12:13,14:15,16:52\n".format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class SomaticTagTestCase(unittest.TestCase):

    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}HC_SOM,Number=1,Type=Integer,Description="Jacquard somatic status for Strelka: 0=non-somatic,1=somatic (based on PASS in FILTER column)",Source="Jacquard",Version={1}>'.format(strelka.JQ_STRELKA_TAG, __version__), strelka._SomaticTag().metaheader)

    def test_format_missingPASS(self):
        tag = strelka._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:{0}HC_SOM|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_presentPASS(self):
        tag = strelka._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|PASS|INFO|SS:F2:F3|2:SA.2:SA.3|5:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|PASS|INFO|SS:F2:F3:{0}HC_SOM|2:SA.2:SA.3:0|5:SB.2:SB.3:1\n").format(strelka.JQ_STRELKA_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA","SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())


class StrelkaTestCase(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        self.caller = strelka.Strelka()

    def test_validate_vcfs_in_directory(self):
        in_files = ["A.vcf","B.vcf"]
        self.caller.validate_vcfs_in_directory(in_files)

        in_files = ["A.vcf","B"]
        self.assertRaisesRegexp(JQException, "ERROR: Non-VCF file in directory. Check parameters and try again", self.caller.validate_vcfs_in_directory, in_files)

    def test_decorate_files(self):
        filenames = ["A/A.strelka.snvs.vcf","A.strelka.indels.vcf"]
        decorator = "normalized"

        actual_filenames = self.caller.decorate_files(filenames, decorator)
        expected_filenames = "A.strelka.normalized.vcf"

        self.assertEquals(expected_filenames,actual_filenames)

        filenames = ["A.strelka.vcf","A.strelka.vcf"]
        decorator = "normalized"

        self.assertRaisesRegexp(JQException, "Each patient in a Strelka directory should have a snvs file and an indels file.", self.caller.decorate_files, filenames, decorator)

    def test_normalize(self):
        writer = MockWriter()
        content1 = ["##foo", "##bar", "#baz"]
        content2 = ["##hi", "##bar", "#baz"]
        reader1 = MockFileReader("indels.vcf", content1)
        reader2 = MockFileReader("snvs.vcf", content2)
        self.caller.normalize(writer,[reader1,reader2])

        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)
        self.assertEquals(["##bar", "##foo", "##hi", "#baz"], writer.lines())

    def test_normalize_mismatchedColumnHeaders(self):
        writer = MockWriter()
        content1 = ["##foo", "##bar", "#baz"]
        content2 = ["##foo", "##bar", "#bluh"]
        reader1 = MockFileReader("indels.vcf", content1)
        reader2 = MockFileReader("snvs.vcf", content2)

        self.assertRaisesRegexp(JQException, r"The column headers for VCF files \[indels.vcf,snvs.vcf\] do not match.",
                                 self.caller.normalize, writer, [reader1,reader2])

    def test_normalize_wrongNumberOfFiles(self):
        self.assertRaisesRegexp(JQException,
                                r"Strelka directories should have exactly two input files per patient, but found \[1\].",
                                self.caller.normalize, MockWriter(), [MockFileReader(input_filepath="foo")])

    def test_normalize_raisesExceptionMissingIndelSnvs(self):
        self.assert_two_vcf_files_throw_exception("foo", "bar")
        self.assert_two_vcf_files_throw_exception("snvs", "bar")
        self.assert_two_vcf_files_throw_exception("foo.snvs", "bar")
        self.assert_two_vcf_files_throw_exception("foo.indels", "bar")
        self.assert_two_vcf_files_throw_exception("foo.indels", "bar.indels")
        self.assert_two_vcf_files_throw_exception("foo.snvs", "bar.snvs")
        self.assert_two_vcf_files_throw_exception("snvs/foo", "indels/bar")
        self.assert_two_vcf_files_throw_exception("indels.snvs", "bar")
        self.assert_two_vcf_files_throw_exception("A.indels.snvs", "B.indels.snvs")

    def test_normalize_writesSequentialRecords(self):
        writer = MockWriter()
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        record2 = "chr2\t.\t.\t.\t.\t.\t.\t.\t."
        record3 = "chr3\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "#bar", record2, record3]
        content2 = ["##foo", "#bar", record1, record3]
        reader1 = MockFileReader("indels.vcf", content1)
        reader2 = MockFileReader("snvs.vcf", content2)
        self.caller.normalize(writer, [reader1, reader2])

        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)
        self.assertEquals(["##foo", "#bar", record1, record2, record3, record3], writer.lines())

    def assert_two_vcf_files_throw_exception(self, file1, file2):
        with self.assertRaisesRegexp(JQException,
                                     r"Each patient in a Strelka directory should have a snvs file and an indels file."):
            self.caller.normalize(MockWriter(), [MockFileReader(input_filepath=file1),
                                                 MockFileReader(input_filepath=file2)])

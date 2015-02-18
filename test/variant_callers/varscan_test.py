# pylint: disable=line-too-long,too-many-public-methods,too-few-public-methods
# pylint: disable=invalid-name,global-statement
import os
import unittest

from testfixtures import TempDirectory

from jacquard import __version__
from jacquard.utils import JQException
from jacquard.variant_callers import varscan
import jacquard.variant_callers.common_tags as common_tags
import jacquard.vcf as vcf
import test.test_case as test_case
from test.vcf_test import MockVcfReader
from jacquard.variant_callers.varscan import _HCTag


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

    def open (self):
        self.opened = True

    def write(self, content):
        self._content.extend(content.splitlines())

    def lines(self):
        return self._content

    def close(self):
        self.closed = True

class MockFileReader(object):
    def __init__(self, input_filepath="/foo/mockFileReader.txt", content = None):
        self.input_filepath = input_filepath
        self.file_name = os.path.basename(input_filepath)
        if content:
            self._content = content
        else:
            self._content = []
        self.open_was_called = False
        self.close_was_called = False

    def open(self):
        self.open_was_called = True

    def read_lines(self):
        for line in self._content:
            yield line

    def close(self):
        self.close_was_called = True

class CommonTagTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(CommonTagTestCase, self).setUp()
        global ORIGINAL_REPORTED_TAG
        global ORIGINAL_PASSED_TAG
        ORIGINAL_REPORTED_TAG = common_tags.ReportedTag
        ORIGINAL_PASSED_TAG = common_tags.PassedTag
        common_tags.ReportedTag = MockCommonTag
        common_tags.PassedTag = MockCommonTag

    def tearDown(self):
        common_tags.ReportedTag = ORIGINAL_REPORTED_TAG
        common_tags.PassedTag = ORIGINAL_PASSED_TAG
        super(CommonTagTestCase, self).tearDown()

    def test_reported_tag(self):
        varscan_instance = varscan.Varscan()
        reported_tag = varscan_instance.tags[0]
        passed_tag = varscan_instance.tags[1]
        self.assertEquals("JQ_VS_", reported_tag.input_caller_name)
        self.assertEquals("JQ_VS_", passed_tag.input_caller_name)

class HCTagTestCase(test_case.JacquardBaseTestCase):
    def test_add_tag_values_highConfidence(self):
        record = vcf.VcfRecord("chr1", "42", "ref", "alt", vcf_filter="pass")
        input_reader = MockFileReader("foo.txt", ["chrom\tposition\tref\tvar","chr1\t42\tref\tvar","chr2\t50\tref\tvar"])
        expected = "pass"
        actual = _HCTag(input_reader).add_tag_values(record).filter
        self.assertEquals(expected, actual)

    def test_add_tag_values_lowConfidence(self):
        record = vcf.VcfRecord("chr1", "30", "ref", "alt", vcf_filter="pass")
        input_reader = MockFileReader("foo.txt", ["chrom\tposition\tref\tvar","chr1\t42\tref\tvar","chr2\t50\tref\tvar"])
        expected = varscan._LOW_CONFIDENCE_FILTER
        actual = _HCTag(input_reader).add_tag_values(record).filter
        self.assertEquals(expected, actual)

    def test_add_tag_values_appendsLowConfidenceForNonPass(self):
        record = vcf.VcfRecord("chr1", "30", "ref", "alt", vcf_filter="indelBias")
        input_reader = MockFileReader("foo.txt", ["chrom\tposition\tref\tvar","chr1\t42\tref\tvar","chr2\t50\tref\tvar"])
        expected = "indelBias;" + varscan._LOW_CONFIDENCE_FILTER
        actual = _HCTag(input_reader).add_tag_values(record).filter
        self.assertEquals(expected, actual)

class AlleleFreqTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}AF,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ)",Source="Jacquard",Version={1}>'.format(varscan.JQ_VARSCAN_TAG, __version__),
                         varscan._AlleleFreqTag().metaheader)

    def test_format_missingAFTag(self):
        tag = varscan._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        originalVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())

    def test_format_presentAFTag(self):
        tag = varscan._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3|56.7%:SA.2:SA.3|83.4%:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3:{0}AF|56.7%:SA.2:SA.3:0.57|83.4%:SB.2:SB.3:0.83\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_multAlt(self):
        tag = varscan._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3|56.7%,41%:SA.2:SA.3|83.4%,23%:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3:{0}AF|56.7%,41%:SA.2:SA.3:0.57,0.41|83.4%,23%:SB.2:SB.3:0.83,0.23\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}DP,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP)",Source="Jacquard",Version={1}>'.format(varscan.JQ_VARSCAN_TAG, __version__),
                         varscan._DepthTag().metaheader)

    def test_format_missingDPTag(self):
        tag = varscan._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        originalVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())

    def test_format_presentDPTag(self):
        tag = varscan._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3|2:SA.2:SA.3|4:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3:{0}DP|2:SA.2:SA.3:2|4:SB.2:SB.3:4\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}HC_SOM,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1=somatic (based on SOMATIC info tag and if sample is TUMOR)",Source="Jacquard",Version={1}>'.format(varscan.JQ_VARSCAN_TAG, __version__),
                         varscan._SomaticTag().metaheader)

    def test_format_missingSSTag(self):
        tag = varscan._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:{0}HC_SOM|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA","SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_presentSSTag_withHC(self):
        tag = varscan._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|{0}HC;SS=2|F1:F2:F3|2:SA.2:SA.3|SB.1:SB.2:SB.3\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|{0}HC;SS=2|F1:F2:F3:{0}HC_SOM|2:SA.2:SA.3:0|SB.1:SB.2:SB.3:1\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_presentSSTag_withoutHC(self):
        tag = varscan._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|SS=2|F1:F2:F3|2:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|SS=2|F1:F2:F3:{0}HC_SOM|2:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_SSTag_notEqual2(self):
        tag = varscan._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|SS=5|F1:F2:F3|2:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|SS=5|F1:F2:F3:{0}HC_SOM|2:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class VarscanTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(VarscanTestCase, self).setUp()
        self.caller = varscan.Varscan()

    def test_validate_vcfs_in_directory(self):
        in_files = ["A.vcf",
                    "B.vcf",
                    "A.somatic.hc.fpfilter.pass"]
        self.caller.validate_vcfs_in_directory(in_files)

        in_files = ["A.vcf", "B"]
        self.assertRaisesRegexp(JQException, "ERROR: Non-VCF or fpfilter file in directory. Check parameters and try again", self.caller.validate_vcfs_in_directory, in_files)

    def test_decorate_files(self):
        filenames = ["A/A.varscan.snp.vcf",
                     "A.varscan.indel.vcf",
                     "A.snp.somatic.hc.fpfilter.pass",
                     "A.indel.somatic.hc.fpfilter.pass"]
        decorator = "normalized"

        actual_filenames = self.caller.decorate_files(filenames, decorator)
        expected_filenames = "A.varscan.normalized.vcf"

        self.assertEquals(expected_filenames,actual_filenames)

        filenames = ["A.varscan.vcf","A.varscan.vcf"]
        decorator = "normalized"

        self.assertRaisesRegexp(JQException, "Each patient in a VarScan directory should have a snp file and an indel file.", self.caller.decorate_files, filenames, decorator)


    def test_normalize(self):
        writer = MockWriter()
        content1 = ["##foo", "##bar", "#baz"]
        content2 = ["##hi", "##bar", "#baz"]
        readers = []
        readers.append(MockFileReader("indel.vcf", content1))
        readers.append(MockFileReader("snp.vcf", content2))
        self.append_hc_files(readers)
        self.caller.normalize(writer, readers)

        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)
        self.assertEquals(["##bar", "##foo", "##hi", "#baz"], writer.lines())

    def test_normalize_modifyBasedOnHCVariants(self):
        writer = MockWriter()
        column_header = "#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleA".replace("|", "\t")
        content1 = ["##FOO\n", column_header, "chr1|161332554|.|A|G|.|.|.|.|.\n".replace("|", "\t")]
        content2 = ["##FOO\n", column_header, "chr1|161332554|.|G|G|.|.|.|.|.\n".replace("|", "\t")]
        readers = []
        readers.append(MockFileReader("indel.vcf", content1))
        readers.append(MockFileReader("snp.vcf", content2))
        self.append_hc_files(readers, content1=["chr1|161332554|A|G|25|1|3.85%|A|25|9|26.47%|R|Somatic|1.0|0.019827310521266846|12|13|3|6|15|10|0|1\n".replace("|", "\t")])
        self.caller.normalize(writer, readers)
        self.maxDiff = None

        self.assertEquals(["##FOO",
                           '##INFO=<ID={0}HC'
                           ',Number=1,Type=Flag,Description="Jacquard high-confidence '
                           'somatic flag for VarScan. Based on intersection with '
                           'filtered VarScan variants.">'.format(varscan.JQ_VARSCAN_TAG),
                           column_header,
                           "chr1|161332554|.|A|G|.|.|.;{0}HC|.|.".format(varscan.JQ_VARSCAN_TAG).replace("|", "\t"),
                           "chr1|161332554|.|G|G|.|.|.|.|.".replace("|", "\t")], writer.lines())

    def test_normalize_wrongNumberOfFiles(self):
        self.assertRaisesRegexp(JQException,
                                r"VarScan directories should have exactly two input VCF files per patient, but found \[1\].",
                                self.caller.normalize,
                                MockWriter(),
                                [MockFileReader(input_filepath="foo.vcf", content=["##foo", "#bar"])])

    def test_normalize_mismatchedColumnHeaders(self):
        writer = MockWriter()
        content1 = ["##foo", "##bar", "#baz"]
        content2 = ["##foo", "##bar", "#bluh"]
        readers = []
        readers.append(MockFileReader("indel.vcf", content1))
        readers.append(MockFileReader("snp.vcf", content2))
        self.append_hc_files(readers)

        self.assertRaisesRegexp(JQException, r"The column headers for VCF files \[indel.vcf,snp.vcf\] do not match.",
                                 self.caller.normalize, writer, readers)

    def test_process_hc_files(self):
        hc_readers = []
        with TempDirectory() as input_file:
            input_file.write("A",
                            ("chrom|position|ref|var|normal_reads1|"+\
                            "normal_reads2|normal_var_freq|normal_gt|"+\
                            "tumor_reads1|tumor_reads2|tumor_var_freq|"+\
                            "tumor_gt|somatic_status|variant_p_value|"+\
                            "somatic_p_value|tumor_reads1_plus|tumor_reads1_minus|"+\
                            "tumor_reads2_plus|tumor_reads2_minus|normal_reads1_plus|"+\
                            "normal_reads1_minus|normal_reads2_plus|normal_reads2_minus\n"+\
                            "chr1|161332554|A|G|25|1|3.85%|A|25|9|26.47%|R|Somatic|1.0|0.019827310521266846|12|13|3|6|15|10|0|1\n"+\
                            "chr2|161332557|G|A|25|1|3.85%|A|25|9|26.47%|R|Somatic|1.0|0.019827310521266846|12|13|3|6|15|10|0|1\n"+\
                            "chr3|99463179|G|A|22|0|3.85%|A|25|9|26.47%|R|Somatic|1.0|0.019827310521266846|12|13|3|6|15|10|0|1\n"+
                            "").replace("|","\t"))
            hc_readers.append(vcf.FileReader(os.path.join(input_file.path,"A")))

            actual = self.caller._process_hc_files(hc_readers)

            metaheader = ('##INFO=<ID={0}HC'
                          ',Number=1,Type=Flag,Description="Jacquard '
                          'high-confidence somatic flag for VarScan. Based on '
                          'intersection with filtered VarScan '
                          'variants.">').format(varscan.JQ_VARSCAN_TAG)
            hc_records = [vcf.VcfRecord("chr1", "161332554", "A", "G"),
                          vcf.VcfRecord("chr2", "161332557", "G", "A"),
                          vcf.VcfRecord("chr3", "99463179", "G", "A")]
            expected = (metaheader, hc_records)
            self.assertEquals(expected, actual)

    def test_normalize_raisesExceptionMissingIndelSnvs(self):
        self.assert_two_vcf_files_throw_exception("foo.vcf", "bar.vcf")
        self.assert_two_vcf_files_throw_exception("snp.vcf", "bar.vcf")
        self.assert_two_vcf_files_throw_exception("foo.snp.vcf", "bar.vcf")
        self.assert_two_vcf_files_throw_exception("foo.indel.vcf", "bar.vcf")
        self.assert_two_vcf_files_throw_exception("foo.indel.vcf", "bar.indel.vcf")
        self.assert_two_vcf_files_throw_exception("foo.snp.vcf", "bar.snp.vcf")
        self.assert_two_vcf_files_throw_exception("snp/foo.vcf", "indel/bar.vcf")
        self.assert_two_vcf_files_throw_exception("indel.snp.vcf", "bar.vcf")
        self.assert_two_vcf_files_throw_exception("A.indel.snp.vcf", "B.indel.snp.vcf")

    def test_normalize_writesSequentialRecords(self):
        writer = MockWriter()
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        record2 = "chr2\t.\t.\t.\t.\t.\t.\t.\t."
        record3 = "chr3\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "#bar", record2, record3]
        content2 = ["##foo", "#bar", record1, record3]
        readers = []
        readers.append(MockFileReader("indel.vcf", content1))
        readers.append(MockFileReader("snp.vcf", content2))
        self.append_hc_files(readers)
        self.caller.normalize(writer,readers)

        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)
        self.assertEquals(["##foo", "#bar", record1, record2, record3, record3], writer.lines())

    def test_normalize_missingHCFilesRaisesException(self):
        writer = MockWriter()
        content1 = ["##foo", "#bar"]
        content2 = ["##foo", "#bar"]
        readers = []
        readers.append(MockFileReader("indel.vcf", content1))
        readers.append(MockFileReader("snp.vcf", content2))
        readers.append(MockFileReader("snp.somatic.hc.fpfilter.pass", []))
        self.assertRaisesRegexp(JQException, r"VarScan directories should have exactly 2 input somatic fpfilter files per patient, but found \[1\].",
                                self.caller.normalize, writer, readers)
        readers = []
        readers.append(MockFileReader("indel.vcf", content1))
        readers.append(MockFileReader("snp.vcf", content2))
        readers.append(MockFileReader("foo", []))
        readers.append(MockFileReader("bar", []))
        self.assertRaisesRegexp(JQException, r"VarScan directories should have exactly 2 input somatic fpfilter files per patient, but found \[0\].",
                                self.caller.normalize, writer, readers)

    def test_normalize_wrongHCFilesRaisesException(self):
        self.assert_two_vcf_files_throw_exception("foo.vcf", "bar.vcf")
        self.assert_two_vcf_files_throw_exception("snp.vcf", "bar.vcf")
        self.assert_two_vcf_files_throw_exception("foo.snp.vcf", "bar.vcf")
        self.assert_two_vcf_files_throw_exception("foo.indel.vcf", "bar.vcf")
        self.assert_two_vcf_files_throw_exception("foo.indel.vcf", "bar.indel.vcf")
        self.assert_two_vcf_files_throw_exception("foo.snp.vcf", "bar.snp.vcf")
        self.assert_two_vcf_files_throw_exception("snp/foo.vcf", "indel/bar.vcf")
        self.assert_two_vcf_files_throw_exception("indel.snp.vcf", "bar.vcf")
        self.assert_two_vcf_files_throw_exception("A.indel.snp.vcf", "B.indel.snp.vcf")

    def assert_two_hc_files_throw_exception(self, file1, file2):
        writer = MockWriter()
        readers = []
        readers.append(MockFileReader("indel.vcf", ["##foo", "#bar"]))
        readers.append(MockFileReader("snp.vcf", ["##foo", "#bar"]))
        self.append_hc_files(readers, file1, file2)
        with self.assertRaisesRegexp(JQException, r"ERROR: Each patient in a VarScan directory should have a somatic fpfilter snp file and indel file."):
            self.caller.normalize(writer, readers)

    def assert_two_vcf_files_throw_exception(self, file1, file2):
        readers = []
        content = ["##foo","#bar"]
        readers.append(MockFileReader(input_filepath=file1, content=content))
        readers.append(MockFileReader(input_filepath=file2, content=content))
        self.append_hc_files(readers)

        with self.assertRaisesRegexp(JQException,r"Each patient in a VarScan directory should have a snp file and an indel file."):
            self.caller.normalize(MockWriter(), readers)

    @staticmethod
    def append_hc_files(readers, file1="snp.somatic.hc.fpfilter.pass", file2="indel.somatic.hc.fpfilter.pass", content1=None, content2=None):
        if not content1:
            content1 = []
        if not content2:
            content2 = []
        readers.append(MockFileReader(file1, content1))
        readers.append(MockFileReader(file2, content2))

    def test_claim_vcf_only(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##MuTect=123", "#chrom", record1]
        content2 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("fileA.vcf", content1)
        reader2 = MockFileReader("fileB.vcf", content2)
        file_readers = [reader1, reader2]

        caller = varscan.Varscan()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(1, len(unrecognized_readers))
        self.assertEquals([reader1], unrecognized_readers)
        self.assertEquals(1, len(vcf_readers))
        self.assertIsInstance(vcf_readers[0], vcf.RecognizedVcfReader)
        self.assertEquals(reader2.file_name, vcf_readers[0].file_name)

    def test_claim_vcf_and_filter_file(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = [self.entab("chrom|position|ref|var"),
                    record1]
        content2 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("patientA.indel.Somatic.hc.fpfilter.pass", content1)
        reader2 = MockFileReader("patientA.indel.vcf", content2)
        reader3 = MockFileReader("patientA.snp.Somatic.hc.fpfilter.pass", content1)
        reader4 = MockFileReader("patientA.snp.vcf", content2)
        reader5 = MockFileReader("patientA.readme", ["foo"])
        file_readers = [reader1, reader2, reader3, reader4, reader5]

        caller = varscan.Varscan()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(1, len(unrecognized_readers))
        self.assertEquals([reader5], unrecognized_readers)
        self.assertEquals(2, len(vcf_readers))
        self.assertIsInstance(vcf_readers[0], vcf.RecognizedVcfReader)
        self.assertEquals(reader2.file_name, vcf_readers[0].file_name)
        self.assertIsInstance(vcf_readers[1], vcf.RecognizedVcfReader)
        self.assertEquals(reader4.file_name, vcf_readers[1].file_name)


    def test_claim_ignores_unpaired_non_vcf_files(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("fileA.txt", content1)
        file_readers = [reader1]

        caller = varscan.Varscan()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(1, len(unrecognized_readers))
        self.assertEquals([reader1], unrecognized_readers)
        self.assertEquals(0, len(vcf_readers))

class VarscanVcfReaderTestCase(test_case.JacquardBaseTestCase):
    def test_metaheaders(self):
        vcf_reader = MockVcfReader(metaheaders=["##foo", "##source=VarScan2"])
        varscan_vcf_reader = varscan._VarscanVcfReader(vcf_reader)
        metaheaders = varscan_vcf_reader.metaheaders

        self.assertIn(varscan._AlleleFreqTag().metaheader, metaheaders)
        self.assertIn(varscan._DepthTag().metaheader, metaheaders)
        self.assertIn(varscan._SomaticTag().metaheader, metaheaders)
        self.assertIn("##foo", metaheaders)
        self.assertIn("##source=VarScan2", metaheaders)

    def test_vcf_records_newTagsPresent(self):
        record1 = vcf.VcfRecord(chrom="chr1",
                                pos="21",
                                ref="A",
                                alt="G",
                                sample_tag_values={"sampleA": {"DP": "45"},
                                                   "sampleB": {"DP": "67"}})
        record2 = vcf.VcfRecord(chrom="chr1",
                                pos="22",
                                ref="A",
                                alt="G",
                                sample_tag_values={"sampleA": {"DP": "46"},
                                                   "sampleB": {"DP": "68"}})
        vcf_reader = MockVcfReader(records=[record1, record2])

        varscan_vcf_reader = varscan._VarscanVcfReader(vcf_reader)
        vcf_records = [record for record in varscan_vcf_reader.vcf_records()]

        self.assertEquals(2, len(vcf_records))
        self.assertIn(varscan.JQ_VARSCAN_TAG + "DP",
                      vcf_records[0].sample_tag_values["sampleA"])
        self.assertIn(varscan.JQ_VARSCAN_TAG + "DP",
                      vcf_records[1].sample_tag_values["sampleA"])

#TODO: jebene - this isn't an accurate reflection of a match - MAKE PASS
    def xtest_vcf_records_SomHcFileIndel(self):
        record1 = vcf.VcfRecord(chrom="chr1", pos="21", ref="A", alt="G", vcf_filter="PASS")
        record2 = vcf.VcfRecord(chrom="chr1", pos="22", ref="A", alt="TT", vcf_filter="PASS")
        vcf_reader = MockVcfReader(records=[record1, record2])

        content1 = ["chrom\tposition\tref\tvar",
                    "chr1\t21\tT\t-C",
                    "chr1\t22\tA\t+TT"]
        somatic_hc_reader = MockFileReader("fileA.Somatic.hc.fpfilter.pass", content1)

        varscan_vcf_reader = varscan._VarscanVcfReader(vcf_reader, somatic_hc_reader)
        vcf_records = [record for record in varscan_vcf_reader.vcf_records()]

        self.assertEquals(2, len(vcf_records))
        self.assertIn(varscan._LOW_CONFIDENCE_FILTER, vcf_records[0].filter)
        self.assertIn("PASS", vcf_records[1].filter)

    def test_vcf_records_SomHcFileSNP(self):
        record1 = vcf.VcfRecord(chrom="chr1", pos="21", ref="A", alt="G", vcf_filter="PASS")
        record2 = vcf.VcfRecord(chrom="chr1", pos="22", ref="A", alt="T", vcf_filter="PASS")
        vcf_reader = MockVcfReader(records=[record1, record2])

        content1 = ["chrom\tposition\tref\tvar",
                    "chr1\t21\tA\tG",
                    "chr1\t22\tA\tT"]
        somatic_hc_reader = MockFileReader("fileA.Somatic.hc.fpfilter.pass", content1)

        varscan_vcf_reader = varscan._VarscanVcfReader(vcf_reader, somatic_hc_reader)
        vcf_records = [record for record in varscan_vcf_reader.vcf_records()]

        self.assertEquals(2, len(vcf_records))
        self.assertIn("PASS", vcf_records[0].filter)
        self.assertIn("PASS", vcf_records[1].filter)

    def test_open_and_close(self):
        vcf_reader = MockVcfReader(metaheaders=["##foo", "##source=VarScan2"])
        varscan_vcf_reader = varscan._VarscanVcfReader(vcf_reader)
        varscan_vcf_reader.open()
        varscan_vcf_reader.close()

        self.assertTrue(varscan_vcf_reader.open)
        self.assertTrue(varscan_vcf_reader.close)

# pylint: disable=line-too-long,too-many-public-methods,too-few-public-methods
# pylint: disable=invalid-name,global-statement
from __future__ import print_function, absolute_import, division

from argparse import Namespace
from collections import OrderedDict
import re

import jacquard.utils.logger
import jacquard.utils.utils as utils
import jacquard.utils.vcf as vcf
from jacquard.variant_caller_transforms import varscan
from jacquard.variant_caller_transforms.varscan import _HCTag
import test.utils.mock_logger
import test.utils.test_case as test_case
from test.utils.vcf_test import MockFileReader, MockVcfReader, MockVcfRecord

class HCTagTestCase(test_case.JacquardBaseTestCase):
    def test_add_tag_values_highConfidenceDoesNotChangeFilter(self):
        record = vcf.VcfRecord("chr1", "42", "ref", "alt", vcf_filter="pass")
        input_reader = MockFileReader("foo.txt", ["chrom\tposition\tref\tvar",
                                                  "chr1\t42\tref\tvar",
                                                  "chr2\t50\tref\tvar"])
        expected = "pass"
        actual = _HCTag(input_reader).add_tag_values(record).filter
        self.assertEquals(expected, actual)

    def test_add_tag_values_lowConfidencePassingReplacesFilter(self):
        record = vcf.VcfRecord("chr1", "30", "ref", "alt", vcf_filter="pass")
        input_reader = MockFileReader("foo.txt", ["chrom\tposition\tref\tvar",
                                                  "chr1\t42\tref\tvar",
                                                  "chr2\t50\tref\tvar"])
        tag = _HCTag(input_reader)
        expected = tag._TAG_ID
        actual = tag.add_tag_values(record).filter
        self.assertEquals(expected, actual)


class GenotypeTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}GT,Number=1,Type=String,Description="Jacquard genotype (based on GT)">'.format(varscan.JQ_VARSCAN_TAG),
                         varscan._GenotypeTag().metaheader)

    def test_add_tag_values_addsGTTag(self):
        vcf_record = MockVcfRecord("chr1", "12", "A", "T", vcf_format="AF:GT", samples=["0.2:0/0", "0.4:0/1"])
        tag = varscan._GenotypeTag()
        tag.add_tag_values(vcf_record)

        expected_sample1 = OrderedDict(sorted({"AF": "0.2",
                                               "GT": "0/0",
                                               "{}GT".format(varscan.JQ_VARSCAN_TAG): "0/0"}.items()))
        self.assertEquals(expected_sample1, vcf_record.sample_tag_values[0])

        expected_sample2 = OrderedDict(sorted({"AF": "0.4",
                                               "GT": "0/1",
                                               "{}GT".format(varscan.JQ_VARSCAN_TAG): "0/1"}.items()))
        self.assertEquals(expected_sample2, vcf_record.sample_tag_values[1])

    def test_add_tag_values_missingGTTag(self):
        vcf_record = MockVcfRecord("chr1", "12", "A", "T", vcf_format="AF", samples=["0.2", "0.4"])
        tag = varscan._GenotypeTag()
        tag.add_tag_values(vcf_record)

        expected_sample1 = OrderedDict(sorted({"AF": "0.2"}.items()))
        self.assertEquals(expected_sample1, vcf_record.sample_tag_values[0])

        expected_sample2 = OrderedDict(sorted({"AF": "0.4"}.items()))
        self.assertEquals(expected_sample2, vcf_record.sample_tag_values[1])

class AlleleFreqTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}AF,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ)">'.format(varscan.JQ_VARSCAN_TAG),
                         varscan._AlleleFreqTag().metaheader)

    def test_format_missingAFTag(self):
        tag = varscan._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        originalVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.text(), processedVcfRecord.text())

    def test_format_presentAFTag(self):
        tag = varscan._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3|56.7%:SA.2:SA.3|83.4%:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3:{0}AF|56.7%:SA.2:SA.3:0.57|83.4%:SB.2:SB.3:0.83\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_multAlt(self):
        tag = varscan._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3|56.7%,41%:SA.2:SA.3|83.4%,23%:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FREQ:F2:F3:{0}AF|56.7%,41%:SA.2:SA.3:0.57,0.41|83.4%,23%:SB.2:SB.3:0.83,0.23\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())


class DepthTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}DP,Number=1,Type=Integer,Description="Jacquard depth for VarScan (based on DP)">'.format(varscan.JQ_VARSCAN_TAG),
                         varscan._DepthTag().metaheader)

    def test_format_missingDPTag(self):
        tag = varscan._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        originalVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.text(), processedVcfRecord.text())

    def test_format_presentDPTag(self):
        tag = varscan._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3|2:SA.2:SA.3|4:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3:{0}DP|2:SA.2:SA.3:2|4:SB.2:SB.3:4\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())


class SomaticTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}HC_SOM,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1=somatic (based on SOMATIC info tag and if sample is TUMOR)">'.format(varscan.JQ_VARSCAN_TAG),
                         varscan._SomaticTag().metaheader)

    def test_format_missingSSTag(self):
        tag = varscan._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:{0}HC_SOM|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_presentSSTag_withHC(self):
        tag = varscan._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|PASS|SS=2|F1:F2:F3|2:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|PASS|SS=2|F1:F2:F3:{0}HC_SOM|2:SA.2:SA.3:0|SB.1:SB.2:SB.3:1\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_presentSSTag_withoutHCManyFilters(self):
        tag = varscan._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|{0};foo|SS=2|F1:F2:F3|2:SA.2:SA.3|SB.1:SB.2:SB.3\n".format(varscan._HCTag._TAG_ID).replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|{0};foo|SS=2|F1:F2:F3:{1}HC_SOM|2:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n".format(varscan._HCTag._TAG_ID, varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_presentSSTag_withoutHC(self):
        tag = varscan._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|{0}|SS=2|F1:F2:F3|2:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t").format(varscan._HCTag._TAG_ID)
        expected = "CHROM|POS|ID|REF|ALT|QUAL|{0}|SS=2|F1:F2:F3:{1}HC_SOM|2:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n".format(varscan._HCTag._TAG_ID, varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_SSTag_notEqual2(self):
        tag = varscan._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|PASS|SS=5|F1:F2:F3|2:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|PASS|SS=5|F1:F2:F3:{0}HC_SOM|2:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n".format(varscan.JQ_VARSCAN_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())


class VarscanTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(VarscanTestCase, self).setUp()
        self.caller = varscan.Varscan()
        varscan.logger = test.utils.mock_logger

    def tearDown(self):
        test.utils.mock_logger.reset()
        varscan.logger = jacquard.utils.logger
        super(VarscanTestCase, self).tearDown()

    @staticmethod
    def append_hc_files(readers,
                        file1="snp.somatic.hc.fpfilter.pass",
                        file2="indel.somatic.hc.fpfilter.pass",
                        content1=None,
                        content2=None):
        if not content1:
            content1 = []
        if not content2:
            content2 = []
        readers.append(MockFileReader(file1, content1))
        readers.append(MockFileReader(file2, content2))

    @staticmethod
    def _get_tag_class_names(vcf_reader):
        return [tag.__class__.__name__ for tag in vcf_reader.tags]

    def test_validate_vcf_hc_pairs(self):
        self.caller._validate_vcf_hc_pairs([(MockFileReader("A.vcf"), MockFileReader("A.hc")),
                                            (MockFileReader("B.vcf"), MockFileReader("B.hc"))])
        self.ok()

    def test_validate_vcf_hc_pairs_okIfNoHcFiles(self):
        self.caller._validate_vcf_hc_pairs([(MockFileReader("A.vcf"), None),
                                            (MockFileReader("B.vcf"), None)])
        self.ok()

    def test_validate_vcf_hc_pairs_raisesIfMissingHcFiles(self):
        vcf_hc_pairs = [(MockFileReader("A.vcf"), MockFileReader("A.hc")),
                        (MockFileReader("B.vcf"), None)]
        self.assertRaisesRegexp(utils.UsageError,
                                "The VarScan VCF file \[B.vcf\] has no matching high-confidence file.", 
                                self.caller._validate_vcf_hc_pairs,
                                vcf_hc_pairs)

    def test_validate_vcf_hc_pairs_raisesIfNoMissingVcfFiles(self):
        vcf_hc_pairs = [(MockFileReader("A.vcf"), MockFileReader("A.hc")),
                        (None, MockFileReader("B.hc"))]
        self.assertRaisesRegexp(utils.UsageError,
                                "The VarScan high-confidence file \[B.hc\] has no matching VCF file.", 
                                self.caller._validate_vcf_hc_pairs,
                                vcf_hc_pairs)


    def test_get_hc_file_pattern(self):
        args = Namespace(varscan_hc_filter_file_regex="foo.*")
        compiled_regex = varscan.Varscan._get_hc_file_pattern(args)
        self.assertEquals("foo.*", compiled_regex.pattern)

    def test_get_hc_file_pattern_invalidRegex(self):
        args = Namespace(varscan_hc_filter_file_regex="*foo")
        self.assertRaisesRegexp(utils.UsageError,
                                r"The specified regex \[\*foo\] could not be compiled. Review inputs and try again",
                                varscan.Varscan._get_hc_file_pattern,
                                args)

    def test_validate_filter_file_validFile(self):
        file_reader = MockFileReader("p1.hc.fpfilter.pass", ["chrom\tposition"])
        caller = varscan.Varscan()
        valid_reader = caller._validate_filter_file(file_reader)
        self.assertEquals("p1.hc.fpfilter.pass", valid_reader.file_name)

    def test_validate_filter_file_invalidFile(self):
        file_reader = MockFileReader("p1.hc.fpfilter.pass", ["chrom\tpos\tref"])
        caller = varscan.Varscan()
        valid_reader = caller._validate_filter_file(file_reader)
        self.assertEquals(None, valid_reader)

    def test_claim_multiplePatients(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=VarScan2", "#chrom", record1]
        content2 = ["chrom\tposition", "1\t23"]
        reader1 = MockFileReader("p2.fileA.vcf", content1)
        reader2 = MockFileReader("p2.fileA.Somatic.hc.fpfilter.pass", content2)
        reader3 = MockFileReader("p3.fileA.Somatic.hc.fpfilter.pass", content2)
        reader4 = MockFileReader("p3.fileA.vcf", content1)
        file_readers = [reader1, reader2, reader3, reader4]

        caller = varscan.Varscan()
        dummy, vcf_readers = caller.claim(file_readers)

        self.assertEquals(2, len(vcf_readers))
        self.assertIsInstance(vcf_readers[0], varscan._VarscanVcfReader)
        self.assertEquals("p2.fileA.vcf", vcf_readers[0]._vcf_reader.file_name)
        self.assertIn("_HCTag",
                      self._get_tag_class_names(vcf_readers[0]))
        self.assertIn("_HCTag",
                      self._get_tag_class_names(vcf_readers[1]))
        self.assertEquals(reader1.file_name, vcf_readers[0]._vcf_reader.file_name)

    def test_claim_varscanVcfOnly(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("fileA.snp.vcf", content1)
        reader2 = MockFileReader("fileB.snp.vcf", content1)
        file_readers = [reader1, reader2]

        caller = varscan.Varscan()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(0, len(unrecognized_readers))
        self.assertEquals([], unrecognized_readers)
        self.assertEquals(2, len(vcf_readers))
        self.assertIsInstance(vcf_readers[0], varscan._VarscanVcfReader)

        self.assertEquals(reader1.file_name, vcf_readers[0]._vcf_reader.file_name)

    def test_claim_vcfAndFilterFile(self):
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
        self.assertIsInstance(vcf_readers[0], varscan._VarscanVcfReader)
        self.assertEquals(reader2.file_name, vcf_readers[0]._vcf_reader.file_name)
        self.assertEquals(reader1.file_name, vcf_readers[0]._som_hc_file_reader.file_name)
        self.assertIsInstance(vcf_readers[1], varscan._VarscanVcfReader)
        self.assertEquals(reader4.file_name, vcf_readers[1]._vcf_reader.file_name)
        self.assertEquals(reader3.file_name, vcf_readers[1]._som_hc_file_reader.file_name)

    def test_claim_vcfAndFilterFileNameGiven(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = [self.entab("chrom|position|ref|var"),
                    record1]
        content2 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("patientA.indel.Somatic.foo.bar", content1)
        reader2 = MockFileReader("patientA.indel.vcf", content2)
        reader3 = MockFileReader("patientA.snp.Somatic.foo.bar", content1)
        reader4 = MockFileReader("patientA.snp.vcf", content2)
        reader5 = MockFileReader("patientA.readme", ["foo"])
        file_readers = [reader1, reader2, reader3, reader4, reader5]

        caller = varscan.Varscan()
        caller.hc_file_pattern = re.compile("foo.bar$")
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(1, len(unrecognized_readers))
        self.assertEquals([reader5], unrecognized_readers)
        self.assertEquals(2, len(vcf_readers))
        self.assertIsInstance(vcf_readers[0], varscan._VarscanVcfReader)
        self.assertEquals(reader2.file_name, vcf_readers[0]._vcf_reader.file_name)
        self.assertEquals(reader1.file_name, vcf_readers[0]._som_hc_file_reader.file_name)
        self.assertIsInstance(vcf_readers[1], varscan._VarscanVcfReader)
        self.assertEquals(reader4.file_name, vcf_readers[1]._vcf_reader.file_name)
        self.assertEquals(reader3.file_name, vcf_readers[1]._som_hc_file_reader.file_name)

    def test_claim_vcfAndInvalidFilterFile(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = [self.entab("chrom|pos|ref|alt"),
                    record1]
        content2 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("patientA.indel.Somatic.hc.fpfilter.pass", content1)
        reader2 = MockFileReader("patientA.indel.vcf", content2)
        reader3 = MockFileReader("patientA.snp.Somatic.hc.fpfilter.pass", content1)
        reader4 = MockFileReader("patientA.snp.vcf", content2)
        reader5 = MockFileReader("patientA.readme", ["foo"])
        file_readers = [reader1, reader2, reader3, reader4, reader5]

        caller = varscan.Varscan()
        self.assertRaisesRegexp(utils.JQException,
                                r"The \[2\] input files \[.*\] match high-confidence file names, but the file header is invalid or missing. Review inputs and try again.",
                                caller.claim,
                                file_readers)

    def test_claim_vcfAnd6InvalidFilterFiles(self):
        #pylint: disable=too-many-locals
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = [self.entab("chrom|pos|ref|alt"),
                    record1]
        content2 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("patientA.snp.Somatic.hc.fpfilter.pass", content1)
        reader2 = MockFileReader("patientA.snp.vcf", content2)
        reader3 = MockFileReader("patientB.snp.Somatic.hc.fpfilter.pass", content1)
        reader4 = MockFileReader("patientB.snp.vcf", content2)
        reader5 = MockFileReader("patientC.snp.Somatic.hc.fpfilter.pass", content1)
        reader6 = MockFileReader("patientC.snp.vcf", content2)
        reader7 = MockFileReader("patientD.snp.Somatic.hc.fpfilter.pass", content1)
        reader8 = MockFileReader("patientD.snp.vcf", content2)
        reader9 = MockFileReader("patientE.snp.Somatic.hc.fpfilter.pass", content1)
        reader10 = MockFileReader("patientE.snp.vcf", content2)
        reader11 = MockFileReader("patientF.snp.Somatic.hc.fpfilter.pass", content1)
        reader12 = MockFileReader("patientF.snp.vcf", content2)
        file_readers = [reader1, reader2, reader3, reader4, reader5, reader6,
                        reader7, reader8, reader9, reader10, reader11, reader12]

        caller = varscan.Varscan()
        self.assertRaisesRegexp(utils.JQException,
                                r"The \[6\] input files \[.*\(1 file\(s\) omitted\)\] match high-confidence file names, but the file header is invalid or missing. Review inputs and try again.",
                                caller.claim,
                                file_readers)

    def test_claim_filterRegexDoesNotMatch(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("patientA.vcf", content1)
        file_readers = [reader1]

        caller = varscan.Varscan()
        caller.hc_file_pattern = re.compile("foo.bar$")
        self.assertRaisesRegexp(utils.UsageError,
                                r"The VarScan high-confidence filename regex \[foo\.bar\$\] didn't match any files in the input directory. The beginning of the high-confidence filename must exactly match a VCF filename up to the .vcf extension. Review inputs/command options and try again.",
                                caller.claim,
                                file_readers)

    def test_claim_hcHasNoMatchingVCF(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = [self.entab("chrom|position"),
                    record1]
        content2 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("patientA.snp.vcf", content2)
        reader2 = MockFileReader("patientA.snp.Somatic.hc.fpfilter.pass", content1)
        reader3 = MockFileReader("patientB.snp.Somatic.hc.fpfilter.pass", content1)
        file_readers = [reader1, reader2, reader3]

        caller = varscan.Varscan()
        self.assertRaisesRegexp(utils.JQException,
                                r"The VarScan high-confidence file \[patientB.snp.Somatic.hc.fpfilter.pass\] has no matching VCF file.",
                                caller.claim,
                                file_readers)

    def test_claim_VCFHasNoMatchingHc(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = [self.entab("chrom|position"),
                    record1]
        content2 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("patientA.vcf", content2)
        reader2 = MockFileReader("patientB.vcf", content2)
        reader3 = MockFileReader("patientB.snp.Somatic.hc.fpfilter.pass", content1)
        file_readers = [reader1, reader2, reader3]

        caller = varscan.Varscan()
        self.assertRaisesRegexp(utils.UsageError,
                                r"The VarScan VCF file \[patientA.vcf\] has no matching high-confidence file.",
                                caller.claim,
                                file_readers)

    def test_claim_ignoresUnpairedNonVcfFiles(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("fileA.txt", content1)
        file_readers = [reader1]

        caller = varscan.Varscan()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(1, len(unrecognized_readers))
        self.assertEquals([reader1], unrecognized_readers)
        self.assertEquals(0, len(vcf_readers))

    def test_claim_ignoresOtherCallers(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=Foo", "#chrom", record1]
        content2 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("fileA.txt", content1)
        reader2 = MockFileReader("fileA.vcf", content2)
        file_readers = [reader1, reader2]

        caller = varscan.Varscan()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(1, len(unrecognized_readers))
        self.assertEquals([reader1], unrecognized_readers)
        self.assertEquals(1, len(vcf_readers))

    def test_claim_mismatchingSnpIndelFiles(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("fileA.snp.vcf", content1)
        reader2 = MockFileReader("fileA.indel.vcf", content1)
        reader3 = MockFileReader("fileB.indel.vcf", content1)
        file_readers = [reader1, reader2, reader3]

        caller = varscan.Varscan()
        self.assertRaisesRegexp(utils.JQException,
                                r"Some Varscan VCFs were missing either a snp or indel file. Review inputs/command options and try again.",
                                caller.claim,
                                file_readers)
        actual_log_errors = test.utils.mock_logger.messages["ERROR"]
        expected_log_errors = ["VarScan VCF [fileB.indel] has no snp file."]
        self.assertEquals(expected_log_errors, actual_log_errors)

    def test_claim_allSnpOrIndelOkay(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=VarScan2", "#chrom", record1]
        reader1 = MockFileReader("fileA.indel.vcf", content1)
        reader2 = MockFileReader("fileB.indel.vcf", content1)
        file_readers = [reader1, reader2]

        caller = varscan.Varscan()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(0, len(unrecognized_readers))
        self.assertEquals(2, len(vcf_readers))

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
        self.assertIn("##jacquard.translate.caller=VarScan", metaheaders)

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
                                sample_tag_values={"sampleA": {"FREQ": "46%"},
                                                   "sampleB": {"FREQ": "68%"}})
        vcf_reader = MockVcfReader(records=[record1, record2])

        varscan_vcf_reader = varscan._VarscanVcfReader(vcf_reader)
        vcf_records = [record for record in varscan_vcf_reader.vcf_records()]

        self.assertEquals(2, len(vcf_records))

        self.assertIn("DP", vcf_records[0].format_tags)
        self.assertIn(varscan.JQ_VARSCAN_TAG + "DP", vcf_records[0].format_tags)
        self.assertIn(varscan.JQ_VARSCAN_TAG + "HC_SOM", vcf_records[0].format_tags)
        self.assertIn(varscan.JQ_VARSCAN_TAG + "CALLER_REPORTED", vcf_records[0].format_tags)
        self.assertIn(varscan.JQ_VARSCAN_TAG + "CALLER_PASSED", vcf_records[0].format_tags)

        self.assertIn("FREQ", vcf_records[1].format_tags)
        self.assertIn(varscan.JQ_VARSCAN_TAG + "AF", vcf_records[1].format_tags)
        self.assertIn(varscan.JQ_VARSCAN_TAG + "HC_SOM", vcf_records[1].format_tags)
        self.assertIn(varscan.JQ_VARSCAN_TAG + "CALLER_REPORTED", vcf_records[1].format_tags)
        self.assertIn(varscan.JQ_VARSCAN_TAG + "CALLER_PASSED", vcf_records[1].format_tags)


    def test_vcf_records_SomHcFileSNP(self):
        record1 = vcf.VcfRecord(chrom="chr1", pos="21", ref="A", alt="G", vcf_filter="PASS")
        record2 = vcf.VcfRecord(chrom="chr1", pos="22", ref="A", alt="T", vcf_filter="PASS")
        vcf_reader = MockVcfReader(records=[record1, record2])

        content1 = ["chrom\tposition",
                    "chr1\t21",
                    "chr1\t22"]
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

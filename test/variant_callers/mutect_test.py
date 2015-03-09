# pylint: disable=line-too-long,too-many-public-methods,too-few-public-methods
# pylint: disable=invalid-name,global-statement
import jacquard.variant_callers.mutect as mutect
import jacquard.vcf as vcf
import test.test_case as test_case
from test.vcf_test import MockFileReader, MockVcfReader


class AlleleFreqTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}AF,Number=A,Type=Float,Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA)">'.format(mutect.JQ_MUTECT_TAG), mutect._AlleleFreqTag().metaheader)

    def test_format_missingAFTag(self):
        tag = mutect._AlleleFreqTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        originalVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.text(), processedVcfRecord.text())

    def test_format_presentAFTag(self):
        tag = mutect._AlleleFreqTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3|0.567:SA.2:SA.3|0.834:SB.2:SB.3\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3:{0}AF|0.567:SA.2:SA.3:0.57|0.834:SB.2:SB.3:0.83\n".format(mutect.JQ_MUTECT_TAG))
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_multAlt(self):
        tag = mutect._AlleleFreqTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3|0.5,0.8:SA.2:SA.3|0.7,0.6:SB.2:SB.3\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3:{0}AF|0.5,0.8:SA.2:SA.3:0.5,0.8|0.7,0.6:SB.2:SB.3:0.7,0.6\n".format(mutect.JQ_MUTECT_TAG))
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

class DepthTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}DP,Number=1,Type=Float,Description="Jacquard depth for MuTect (based on DP)">'.format(mutect.JQ_MUTECT_TAG), mutect._DepthTag().metaheader)

    def test_format_missingDPTag(self):
        tag = mutect._DepthTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        originalVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.text(), processedVcfRecord.text())

    def test_format_presentDPTag(self):
        tag = mutect._DepthTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3|2:SA.2:SA.3|4:SB.2:SB.3\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3:{0}DP|2:SA.2:SA.3:2|4:SB.2:SB.3:4\n".format(mutect.JQ_MUTECT_TAG))
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

class SomaticTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}HC_SOM,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1=somatic (based on SS FORMAT tag)">'.format(mutect.JQ_MUTECT_TAG), mutect._SomaticTag().metaheader)

    def test_format_missingSSTag(self):
        tag = mutect._SomaticTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:{0}HC_SOM|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").format(mutect.JQ_MUTECT_TAG)
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_presentSSTag(self):
        tag = mutect._SomaticTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|SS:F2:F3|2:SA.2:SA.3|5:SB.2:SB.3\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|SS:F2:F3:{0}HC_SOM|2:SA.2:SA.3:1|5:SB.2:SB.3:0\n").format(mutect.JQ_MUTECT_TAG)
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

class MutectTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(MutectTestCase, self).setUp()
        self.caller = mutect.Mutect()

    def test_validateInputFile_isValid(self):
        metaheaders = ["##MuTect=blah"]
        self.assertTrue(self.caller.validate_input_file(metaheaders, "#column_header"))

    def test_validateInputFile_isNotValid(self):
        metaheaders = ["Foo"]
        self.assertFalse(self.caller.validate_input_file(metaheaders, "#column_header"))

    def test_claim(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=strelka", "#chrom", record1]
        content2 = ["##foo", "##MuTect=123", "#chrom", record1]
        reader1 = MockFileReader("fileA.vcf", content1)
        reader2 = MockFileReader("fileB.vcf", content2)
        file_readers = [reader1, reader2]

        caller = mutect.Mutect()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(1, len(unrecognized_readers))
        self.assertEquals([reader1], unrecognized_readers)
        self.assertEquals(1, len(vcf_readers))
        self.assertIsInstance(vcf_readers[0], mutect._MutectVcfReader)
        self.assertEquals(reader2.file_name, vcf_readers[0]._vcf_reader.file_name)

    def test_claim_ignoresNonVcfExtensions(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##MuTect=123", "#chrom", record1]
        reader1 = MockFileReader("fileA.txt", content1)
        file_readers = [reader1]

        caller = mutect.Mutect()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(1, len(unrecognized_readers))
        self.assertEquals([reader1], unrecognized_readers)
        self.assertEquals(0, len(vcf_readers))

    def test_claim_vcfExtensionCaseInsensitive(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##MuTect=123", "#chrom", record1]
        reader1 = MockFileReader("fileA.VcF", content1)
        file_readers = [reader1]

        caller = mutect.Mutect()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(0, len(unrecognized_readers))
        self.assertEquals(1, len(vcf_readers))

class MutectVcfReaderTestCase(test_case.JacquardBaseTestCase):
    def test_metaheaders(self):
        vcf_reader = MockVcfReader(metaheaders=["##foo", "##MuTect=123"])
        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)
        metaheaders = mutect_vcf_reader.metaheaders

        self.assertIn(mutect._AlleleFreqTag().metaheader, metaheaders)
        self.assertIn(mutect._DepthTag().metaheader, metaheaders)
        self.assertIn(mutect._SomaticTag().metaheader, metaheaders)
        self.assertIn("##foo", metaheaders)
        self.assertIn("##MuTect=123", metaheaders)
        self.assertIn("##jacquard.translate.caller=MuTect", metaheaders)

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
                                sample_tag_values={"sampleA": {"FA": "0.54"},
                                                   "sampleB": {"FA": "0.76"}})
        vcf_reader = MockVcfReader(records=[record1, record2])

        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)
        vcf_records = [record for record in mutect_vcf_reader.vcf_records()]

        self.assertEquals(2, len(vcf_records))

        self.assertIn("DP", vcf_records[0].format_tags)
        self.assertIn(mutect.JQ_MUTECT_TAG + "DP", vcf_records[0].format_tags)
        self.assertIn(mutect.JQ_MUTECT_TAG + "HC_SOM", vcf_records[0].format_tags)
        self.assertIn(mutect.JQ_MUTECT_TAG + "CALLER_REPORTED", vcf_records[0].format_tags)
        self.assertIn(mutect.JQ_MUTECT_TAG + "CALLER_PASSED", vcf_records[0].format_tags)

        self.assertIn("FA", vcf_records[1].format_tags)
        self.assertIn(mutect.JQ_MUTECT_TAG + "AF", vcf_records[1].format_tags)
        self.assertIn(mutect.JQ_MUTECT_TAG + "HC_SOM", vcf_records[1].format_tags)
        self.assertIn(mutect.JQ_MUTECT_TAG + "CALLER_REPORTED", vcf_records[1].format_tags)
        self.assertIn(mutect.JQ_MUTECT_TAG + "CALLER_PASSED", vcf_records[1].format_tags)


    def test_open_and_close(self):
        vcf_reader = MockVcfReader(metaheaders=["##foo", "##MuTect=123"])
        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)
        mutect_vcf_reader.open()
        mutect_vcf_reader.close()

        self.assertTrue(mutect_vcf_reader.open)
        self.assertTrue(mutect_vcf_reader.close)

    def test_column_header_mangleSampleName(self):
        column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|25714|25715")
        meta_header = '##MuTect="123 tumor_sample_name=25715 normal_sample_name=25714"'
        vcf_reader = MockVcfReader(metaheaders=[meta_header],
                                   column_header=column_header)
        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)

        expected_column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR")

        self.assertEquals(expected_column_header, mutect_vcf_reader.column_header)


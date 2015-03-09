# pylint: disable=line-too-long,too-many-public-methods,too-few-public-methods
# pylint: disable=invalid-name,global-statement,too-many-format-args
from jacquard import __version__
from test.vcf_test import MockFileReader, MockVcfReader
import jacquard.variant_callers.strelka as strelka
import jacquard.vcf as vcf
import test.test_case as test_case


class AlleleFreqTagTestCase(test_case.JacquardBaseTestCase):

    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_SK_AF,Number=A,Type=Float,Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth. Uses (TIR tier 2)/DP2 if available, otherwise uses (ACGT tier2 depth) / DP2)",Source="Jacquard",Version=0.3>'.format(strelka.JQ_STRELKA_TAG, __version__), strelka._AlleleFreqTag().metaheader)

    def test_format_missingAFTag(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        originalVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.text(), processedVcfRecord.text())

    def test_format_AUTagWhenMultAlt(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|A,C|QUAL|FILTER|INFO|AU:CU:GU:TU|1,2:3,4:5,6:7,8|9,10:11,12:13,14:15,16\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|A,C|QUAL|FILTER|INFO|AU:CU:GU:TU:{0}AF|1,2:3,4:5,6:7,8:0.1,0.2|9,10:11,12:13,14:15,16:0.19,0.23\n".format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_AFTag_noAlt(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|.|QUAL|FILTER|INFO|AU:CU:GU:TU|1,2:3,4:5,6:7,8|9,10:11,12:13,14:15,16\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|.|QUAL|FILTER|INFO|AU:CU:GU:TU:{0}AF|1,2:3,4:5,6:7,8:.|9,10:11,12:13,14:15,16:.\n".format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_TIRTag(self):
        tag = strelka._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:TIR|10:3,4|20:11,7\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:TIR:{0}AF|10:3,4:0.4|20:11,7:0.35\n".format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())


class DepthTagTestCase(test_case.JacquardBaseTestCase):

    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}DP,Number=1,Type=Float,Description="Jacquard depth for Strelka (uses DP2 if available, otherwise uses ACGT tier2 depth)",Source="Jacquard",Version={1}>'.format(strelka.JQ_STRELKA_TAG, __version__), strelka._DepthTag().metaheader)

    def test_format_missingTag(self):
        tag = strelka._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        originalVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.text(), processedVcfRecord.text())

    def test_format_DP2Tag(self):
        tag = strelka._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:F2:F3|2:SA.2:SA.3|4:SB.2:SB.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP2:F2:F3:{0}DP|2:SA.2:SA.3:2|4:SB.2:SB.3:4\n".format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_AUTag(self):
        tag = strelka._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|AU:CU:GU:TU|1,2:3,4:5,6:7,8|9,10:11,12:13,14:15,16\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|AU:CU:GU:TU:{0}DP|1,2:3,4:5,6:7,8:20|9,10:11,12:13,14:15,16:52\n".format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())


class SomaticTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}HC_SOM,Number=1,Type=Integer,Description="Jacquard somatic status for Strelka: 0=non-somatic,1=somatic (based on PASS in FILTER column)",Source="Jacquard",Version={1}>'.format(strelka.JQ_STRELKA_TAG, __version__), strelka._SomaticTag().metaheader)

    def test_format_missingPASS(self):
        tag = strelka._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|', "\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:{0}HC_SOM|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").format(strelka.JQ_STRELKA_TAG).replace('|', "\t")
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_presentPASS(self):
        tag = strelka._SomaticTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|PASS|INFO|SS:F2:F3|2:SA.2:SA.3|5:SB.2:SB.3\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|PASS|INFO|SS:F2:F3:{0}HC_SOM|2:SA.2:SA.3:0|5:SB.2:SB.3:1\n").format(strelka.JQ_STRELKA_TAG)
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())


class StrelkaTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(StrelkaTestCase, self).setUp()
        self.caller = strelka.Strelka()

    def test_claim(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=strelka", "#chrom", record1]
        content2 = ["##foo", "##MuTect", "#chrom", record1]
        reader1 = MockFileReader("fileA.vcf", content1)
        reader2 = MockFileReader("fileB.vcf", content2)
        file_readers = [reader1, reader2]

        caller = strelka.Strelka()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(1, len(unrecognized_readers))
        self.assertEquals([reader2], unrecognized_readers)
        self.assertEquals(1, len(vcf_readers))
        self.assertIsInstance(vcf_readers[0], strelka._StrelkaVcfReader)
        self.assertEquals(reader1.file_name, vcf_readers[0]._vcf_reader.file_name)

    def test_claim_ignores_non_vcf_files(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=strelka", "#chrom", record1]
        reader1 = MockFileReader("fileA.txt", content1)
        file_readers = [reader1]

        caller = strelka.Strelka()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(1, len(unrecognized_readers))
        self.assertEquals([reader1], unrecognized_readers)
        self.assertEquals(0, len(vcf_readers))


class StrelkaVcfReaderTestCase(test_case.JacquardBaseTestCase):
    def test_metaheaders(self):
        vcf_reader = MockVcfReader(metaheaders=["##foo", "##source=strelka"])
        strelka_vcf_reader = strelka._StrelkaVcfReader(vcf_reader)
        metaheaders = strelka_vcf_reader.metaheaders

        self.assertIn(strelka._AlleleFreqTag().metaheader, metaheaders)
        self.assertIn(strelka._DepthTag().metaheader, metaheaders)
        self.assertIn(strelka._SomaticTag().metaheader, metaheaders)
        self.assertIn("##foo", metaheaders)
        self.assertIn("##source=strelka", metaheaders)
        self.assertIn("##jacquard.translate.caller=Strelka", metaheaders)

    def test_vcf_records_newTagsPresent(self):
        record1 = vcf.VcfRecord(chrom="chr1",
                                pos="21",
                                ref="A",
                                alt="G",
                                sample_tag_values={"sampleA": {"DP2": "45"},
                                                   "sampleB": {"DP2": "67"}})
        record2 = vcf.VcfRecord(chrom="chr1",
                                pos="22",
                                ref="A",
                                alt="G",
                                sample_tag_values={"sampleA": {"TIR": "10,20", "DP2":"100"},
                                                   "sampleB": {"TIR": "15,25", "DP2":"100"}})
        vcf_reader = MockVcfReader(records=[record1, record2])

        strelka_vcf_reader = strelka._StrelkaVcfReader(vcf_reader)
        vcf_records = [record for record in strelka_vcf_reader.vcf_records()]

        self.assertEquals(2, len(vcf_records))

        self.assertIn("DP2", vcf_records[0].format_tags)
        self.assertIn(strelka.JQ_STRELKA_TAG + "DP", vcf_records[0].format_tags)
        self.assertIn(strelka.JQ_STRELKA_TAG + "HC_SOM", vcf_records[0].format_tags)
        self.assertIn(strelka.JQ_STRELKA_TAG + "CALLER_REPORTED", vcf_records[0].format_tags)
        self.assertIn(strelka.JQ_STRELKA_TAG + "CALLER_PASSED", vcf_records[0].format_tags)

        self.assertIn("TIR", vcf_records[1].format_tags)
        self.assertIn(strelka.JQ_STRELKA_TAG + "AF", vcf_records[1].format_tags)
        self.assertIn(strelka.JQ_STRELKA_TAG + "HC_SOM", vcf_records[1].format_tags)
        self.assertIn(strelka.JQ_STRELKA_TAG + "CALLER_REPORTED", vcf_records[1].format_tags)
        self.assertIn(strelka.JQ_STRELKA_TAG + "CALLER_PASSED", vcf_records[1].format_tags)

    def test_open_and_close(self):
        vcf_reader = MockVcfReader(metaheaders=["##foo", "##source=strelka"])
        strelka_vcf_reader = strelka._StrelkaVcfReader(vcf_reader)
        strelka_vcf_reader.open()
        strelka_vcf_reader.close()

        self.assertTrue(strelka_vcf_reader.open)
        self.assertTrue(strelka_vcf_reader.close)

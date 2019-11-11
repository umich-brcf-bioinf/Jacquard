# pylint: disable=line-too-long,too-many-public-methods,too-few-public-methods
# pylint: disable=invalid-name,global-statement
from __future__ import print_function, absolute_import, division

from collections import OrderedDict

import jacquard.variant_caller_transforms.mutect as mutect
import jacquard.utils.utils as utils
import jacquard.utils.vcf as vcf

import test.utils.test_case as test_case
from test.utils.vcf_test import MockFileReader, MockVcfReader, MockVcfRecord


class GenotypeTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}GT,Number=1,Type=String,Description="Jacquard genotype (based on GT)">'.format(mutect.JQ_MUTECT_TAG),
                         mutect._GenotypeTag().metaheader)

    def test_add_tag_values(self):
        vcf_record = MockVcfRecord("chr1", "24", "A", "T", vcf_format="AF:GT", samples=["0.3:0", "0.2:0/1"])
        tag = mutect._GenotypeTag()
        tag.add_tag_values(vcf_record)

        expected_sample1 = OrderedDict(sorted({"AF": "0.3",
                                               "GT": "0",
                                               "{}GT".format(mutect.JQ_MUTECT_TAG): "0/0"}.items()))
        self.assertEquals(expected_sample1, vcf_record.sample_tag_values[0])

        expected_sample2 = OrderedDict(sorted({"AF": "0.2",
                                               "GT": "0/1",
                                               "{}GT".format(mutect.JQ_MUTECT_TAG): "0/1"}.items()))
        self.assertEquals(expected_sample2, vcf_record.sample_tag_values[1])

    def test_add_tag_values_missingGTTag(self):
        vcf_record = MockVcfRecord("chr1", "24", "A", "T", vcf_format="AF", samples=["0.3", "0.2"])
        tag = mutect._GenotypeTag()
        tag.add_tag_values(vcf_record)

        expected_sample1 = OrderedDict(sorted({"AF": "0.3"}.items()))
        self.assertEquals(expected_sample1, vcf_record.sample_tag_values[0])

        expected_sample2 = OrderedDict(sorted({"AF": "0.2"}.items()))
        self.assertEquals(expected_sample2, vcf_record.sample_tag_values[1])

    def test_add_tag_values_noChangeInGT(self):
        vcf_record = MockVcfRecord("chr1", "24", "A", "T", vcf_format="AF:GT", samples=["0.3:0/1", "0.2:0/1"])
        tag = mutect._GenotypeTag()
        tag.add_tag_values(vcf_record)

        expected_sample1 = OrderedDict(sorted({"AF": "0.3",
                                               "GT": "0/1",
                                               "{}GT".format(mutect.JQ_MUTECT_TAG): "0/1"}.items()))
        self.assertEquals(expected_sample1, vcf_record.sample_tag_values[0])

        expected_sample2 = OrderedDict(sorted({"AF": "0.2",
                                               "GT": "0/1",
                                               "{}GT".format(mutect.JQ_MUTECT_TAG): "0/1"}.items()))
        self.assertEquals(expected_sample2, vcf_record.sample_tag_values[1])


class AlleleFreqTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        expected = ('##FORMAT=<ID={0}AF,Number=A,Type=Float,'
                    'Description="Jacquard allele frequency for MuTect: '
                    'Decimal allele frequency rounded to 4 digits')\
                    .format(mutect.JQ_MUTECT_TAG)
        old_style_metaheaders = ['##FORMAT=<ID=FA,...>']
        self.assertEqual(expected + ' (based on FA)">',
                         mutect._AlleleFreqTag(old_style_metaheaders).metaheader)
        new_style_metaheaders = ['##FORMAT=<ID=AF,...>']
        self.assertEqual(expected + ' (based on AF)">',
                         mutect._AlleleFreqTag(new_style_metaheaders).metaheader)

    def test_missingAlleleFreqTagInInputMetaheader(self):
        inconclusive_metaheaders = []
        self.assertRaisesRegexp(utils.JQException,
                               r'could not determine the correct allele frequency FORMAT tag',
                               mutect._AlleleFreqTag,
                               inconclusive_metaheaders)

    def test_format_missingAFTag(self):
        new_style_metaheaders = ['##FORMAT=<ID=AF,...>']
        tag = mutect._AlleleFreqTag(new_style_metaheaders)
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        originalVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(originalVcfRecord.text(), processedVcfRecord.text())

    def test_format_presentAFTag(self):
        new_style_metaheaders = ['##FORMAT=<ID=AF,...>']
        tag = mutect._AlleleFreqTag(new_style_metaheaders)
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|AF:F2:F3|0.34567:SA.2:SA.3|0.76543:SB.2:SB.3\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|AF:F2:F3:{0}AF|0.34567:SA.2:SA.3:0.3457|0.76543:SB.2:SB.3:0.7654\n".format(mutect.JQ_MUTECT_TAG))
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_multAlt(self):
        new_style_metaheaders = ['##FORMAT=<ID=AF,...>']
        tag = mutect._AlleleFreqTag(new_style_metaheaders)
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|AF:F2:F3|0.5,0.8:SA.2:SA.3|0.7,0.6:SB.2:SB.3\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|AF:F2:F3:{0}AF|0.5,0.8:SA.2:SA.3:0.5,0.8|0.7,0.6:SB.2:SB.3:0.7,0.6\n".format(mutect.JQ_MUTECT_TAG))
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

class DepthTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}DP,Number=1,Type=Integer,Description="Jacquard depth for MuTect (based on DP)">'.format(mutect.JQ_MUTECT_TAG), mutect._DepthTag().metaheader)

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

class SomaticTagSSTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}HC_SOM,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1=somatic (based on SS FORMAT tag)">'.format(mutect.JQ_MUTECT_TAG), mutect._SomaticTagSS().metaheader)

    def test_format_missingSSTag(self):
        tag = mutect._SomaticTagSS()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:{0}HC_SOM|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").format(mutect.JQ_MUTECT_TAG)
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_format_presentSSTag(self):
        tag = mutect._SomaticTagSS()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|SS:F2:F3|2:SA.2:SA.3|5:SB.2:SB.3\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|SS:F2:F3:{0}HC_SOM|2:SA.2:SA.3:1|5:SB.2:SB.3:0\n").format(mutect.JQ_MUTECT_TAG)
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

class SomaticTagFilterMutectCallsTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEqual(\
'''##FORMAT=<ID={0}HC_SOM,Number=1,Type=Integer,Description="Jacquard somatic
 status for MuTect: 0=non-somatic,1=somatic (based on FilterMutectCalls setting
 filter to PASS)">'''.replace('\n', '').format(mutect.JQ_MUTECT_TAG),
            mutect._SomaticTagFilterMutectCalls().metaheader)

    def test_filterFailNotSomatic(self):
        tag = mutect._SomaticTagFilterMutectCalls()
        line = self.entab(\
'CHROM|POS|ID|REF|ALT|QUAL|filter_failed|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n')
        expected = self.entab(\
'CHROM|POS|ID|REF|ALT|QUAL|filter_failed|INFO|F1:F2:F3:{0}HC_SOM|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n'\
            ).format(mutect.JQ_MUTECT_TAG)
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_filterNullNotSomatic(self):
        tag = mutect._SomaticTagFilterMutectCalls()
        line = self.entab(\
'CHROM|POS|ID|REF|ALT|QUAL|.|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n')
        expected = self.entab(\
'CHROM|POS|ID|REF|ALT|QUAL|.|INFO|F1:F2:F3:{0}HC_SOM|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n'\
            ).format(mutect.JQ_MUTECT_TAG)
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_filterPassMissingGenotype(self):
        tag = mutect._SomaticTagFilterMutectCalls()
        line = self.entab(\
'chrQ|42|ID|A|G|QUAL|PASS|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n')
        expected = self.entab(\
'chrQ|42|ID|A|G|QUAL|PASS|INFO|F1:F2:F3:{0}HC_SOM|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n'\
            ).format(mutect.JQ_MUTECT_TAG)
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        self.assertRaisesRegexp(utils.JQException,
                               r'Cannot assign somatic status using FilterMutectCalls.*chrQ:42:A:G',
                               tag.add_tag_values,
                               processedVcfRecord)

    def test_filterPassBothSamplesVariants(self):
        tag = mutect._SomaticTagFilterMutectCalls()
        line = self.entab(\
'CHROM|POS|ID|REF|ALT|QUAL|PASS|INFO|F1:F2:F3:GT|SA.1:SA.2:SA.3:0/1|SB.1:SB.2:SB.3:0/2\n')
        expected = self.entab(\
'CHROM|POS|ID|REF|ALT|QUAL|PASS|INFO|F1:F2:F3:GT:{0}HC_SOM|SA.1:SA.2:SA.3:0/1:1|SB.1:SB.2:SB.3:0/2:1\n'\
            ).format(mutect.JQ_MUTECT_TAG)
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())

    def test_filterPassSomeSamplesVariants(self):
        tag = mutect._SomaticTagFilterMutectCalls()
        line = self.entab(\
'CHROM|POS|ID|REF|ALT|QUAL|PASS|INFO|F1:F2:F3:GT|SA.1:SA.2:SA.3:0/0|SB.1:SB.2:SB.3:0/1\n')
        expected = self.entab(\
'CHROM|POS|ID|REF|ALT|QUAL|PASS|INFO|F1:F2:F3:GT:{0}HC_SOM|SA.1:SA.2:SA.3:0/0:0|SB.1:SB.2:SB.3:0/1:1\n'\
            ).format(mutect.JQ_MUTECT_TAG)
        processedVcfRecord = vcf.VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.text())


class MutectTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(MutectTestCase, self).setUp()
        self.caller = mutect.Mutect()

    def test_claim(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo", "##source=strelka", "#chrom", record1]
        content2 = ["##foo", "##MuTect=123", "##FORMAT=<ID=AF,...>", "#chrom", record1]
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
        content1 = ["##foo", "##MuTect=123", "##FORMAT=<ID=AF,...>", "#chrom", record1]
        reader1 = MockFileReader("fileA.VcF", content1)
        file_readers = [reader1]

        caller = mutect.Mutect()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(0, len(unrecognized_readers))
        self.assertEquals(1, len(vcf_readers))

    def test_claim_metaheaderRecognizesMutectV2x(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo",
                    "##MuTect=2.1",
                    "##FORMAT=<ID=FA,...>",
                    "#chrom",
                    record1]
        reader1 = MockFileReader("fileA.vcf", content1)
        file_readers = [reader1]

        caller = mutect.Mutect()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(0, len(unrecognized_readers))
        self.assertEquals(1, len(vcf_readers))

    def test_claim_metaheaderRecognizesMutectV3x(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo",
                    '##GATKCommandLine.MuTect2=<ID=MuTect2,CommandLineOptions="MuTect2 ...">',
                    "##FORMAT=<ID=AF,...>",
                    "#chrom",
                    record1]
        reader1 = MockFileReader("fileB.vcf", content1)
        file_readers = [reader1]

        caller = mutect.Mutect()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(0, len(unrecognized_readers))
        self.assertEquals(1, len(vcf_readers))

    def test_claim_metaheaderRecognizesMutectV4x(self):
        record1 = "chr1\t.\t.\t.\t.\t.\t.\t.\t."
        content1 = ["##foo",
                    '##GATKCommandLine=<ID=Mutect2,CommandLine="Mutect2 ...">',
                    "##FORMAT=<ID=AF,...>",
                    "#chrom",
                    record1]
        reader1 = MockFileReader("fileB.vcf", content1)
        file_readers = [reader1]

        caller = mutect.Mutect()
        unrecognized_readers, vcf_readers = caller.claim(file_readers)

        self.assertEquals(0, len(unrecognized_readers))
        self.assertEquals(1, len(vcf_readers))

class MutectVcfReaderTestCase(test_case.JacquardBaseTestCase):
    def test_common_metaheaders(self):
        vcf_reader = MockVcfReader(metaheaders=["##foo",
                                                "##MuTect=123",
                                                '##FORMAT=<ID=FA,...>'])
        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)
        metaheaders = mutect_vcf_reader.metaheaders

        self.assertIn(mutect._AlleleFreqTag(vcf_reader.metaheaders).metaheader, metaheaders)
        self.assertIn(mutect._DepthTag().metaheader, metaheaders)
        self.assertIn("##foo", metaheaders)
        self.assertIn("##MuTect=123", metaheaders)
        self.assertIn("##jacquard.translate.caller=MuTect", metaheaders)

    def test_SomaticTagSS_metaheaders(self):
        vcf_reader = MockVcfReader(metaheaders=["##foo",
                                                "##MuTect=123",
                                                '##FORMAT=<ID=FA,...>'])
        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)
        metaheaders = mutect_vcf_reader.metaheaders

        self.assertIn(mutect._SomaticTagSS().metaheader, metaheaders)

    def test_SomaticTagFilterMutectCalls_metaheaders(self):
        vcf_reader = MockVcfReader(metaheaders=["##foo",
                                                "##MuTect=123",
                                                '##FORMAT=<ID=FA,...>',
                                                "##source=FilterMutectCalls"])
        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)
        metaheaders = mutect_vcf_reader.metaheaders

        self.assertIn(mutect._SomaticTagFilterMutectCalls().metaheader, metaheaders)


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
        vcf_reader = MockVcfReader(metaheaders=['##FORMAT=<ID=FA,...>'],
                                   records=[record1, record2])

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
        vcf_reader = MockVcfReader(metaheaders=["##foo",
                                                "##MuTect=123",
                                                '##FORMAT=<ID=FA,...>'])
        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)
        mutect_vcf_reader.open()
        mutect_vcf_reader.close()

        self.assertTrue(mutect_vcf_reader.open)
        self.assertTrue(mutect_vcf_reader.close)

    def test_column_header_mangleSampleNameMutect1(self):
        column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|25714|25715")
        metaheaders = ['##MuTect="123 tumor_sample_name=25715 normal_sample_name=25714"',
                       '##FORMAT=<ID=FA,...>',]
        vcf_reader = MockVcfReader(metaheaders=metaheaders,
                                   column_header=column_header)

        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)

        expected_column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR")

        self.assertEquals(expected_column_header, mutect_vcf_reader.column_header)

    def test_column_header_mangleSampleNameMutect2UsesSampleMetalinesIfAvailable(self):
        column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|25714|25715")
        meta_header = '''
##GATKCommandLine=<ID=Mutect2,CommandLine="Mutect2  --tumor-sample A --normal-sample B",Date="recent">'
##FORMAT=<ID=FA,...>
##foo=42
##SAMPLE=<ID=NORMAL,SampleName=25714,File=foo.bam>
##SAMPLE=<ID=TUMOR,SampleName=25715,File=bar.bam>
##baz=42
'''
        vcf_reader = MockVcfReader(metaheaders=meta_header.strip().split('\n'),
                                   column_header=column_header)
        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)

        expected_column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR")

        self.assertEquals(expected_column_header, mutect_vcf_reader.column_header)

    def test_column_header_mangleSampleNameMutect2UsesCommandLineIfNoSampleMetalines(self):
        column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|25714|25715")
        metaheaders = ['##GATKCommandLine=<ID=Mutect2,CommandLine="Mutect2  --tumor-sample 25715 --normal-sample 25714",Date="recent">',
                       '##FORMAT=<ID=FA,...>',]
        vcf_reader = MockVcfReader(metaheaders=metaheaders,
                                   column_header=column_header)
        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)

        expected_column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR")

        self.assertEquals(expected_column_header, mutect_vcf_reader.column_header)


    def test_column_header_mangleSampleNameMutect2IgnoresHelpFlag(self):
        column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|25714|25715")
        metaheaders = ['##GATKCommandLine=<ID=Mutect2,CommandLine="Mutect2  --tumor-sample 25715 --normal-sample 25714 --help false",Date="recent">',
                       '##FORMAT=<ID=FA,...>',]
        vcf_reader = MockVcfReader(metaheaders=metaheaders,
                                   column_header=column_header)

        mutect_vcf_reader = mutect._MutectVcfReader(vcf_reader)

        expected_column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR")

        self.assertEquals(expected_column_header, mutect_vcf_reader.column_header)

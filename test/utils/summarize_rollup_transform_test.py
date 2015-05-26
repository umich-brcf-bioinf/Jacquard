#pylint: disable=too-few-public-methods, invalid-name, line-too-long
#pylint: disable=too-many-instance-attributes, too-many-public-methods
from __future__ import print_function, absolute_import, division

import re

from jacquard.utils.utils import JQException
import jacquard.utils.utils as utils
import jacquard.variant_caller_transforms.common_tags as common_tags
import jacquard.variant_caller_transforms.mutect as mutect
import jacquard.utils.summarize_rollup_transform as summarize_caller
import jacquard.variant_caller_transforms.varscan as varscan
from jacquard.utils.vcf import VcfRecord
import test.utils.test_case as test_case


class CallersReportedListTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = summarize_caller._CallersReportedListTag().metaheader.split("\n")
        self.assertEquals('##FORMAT=<ID={}CALLERS_REPORTED_LIST,Number=.,Type=String,Description="Comma-separated list variant callers which listed this variant in the Jacquard tagged VCF">'.format(summarize_caller.JQ_SUMMARY_TAG),
                          split_metaheader[0])

    def test_add_tag_values(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:1:1|Y:1:1\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._CallersReportedListTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}CALLERS_REPORTED_LIST|X:1:1:MT,VS|Y:1:1:MT,VS\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED, summarize_caller.JQ_SUMMARY_TAG))
        self.assertEquals(expected, processedVcfRecord.text())

    def test_add_tag_values_nullValues(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:1:.|Y:1:.\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._CallersReportedListTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}CALLERS_REPORTED_LIST|X:1:.:MT|Y:1:.:MT\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED, summarize_caller.JQ_SUMMARY_TAG))
        self.assertEquals(expected, processedVcfRecord.text())

class CallersReportedTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = summarize_caller._CallersReportedTag().metaheader.split("\n")
        self.assertEquals('##FORMAT=<ID={}{},Number=1,Type=Integer,Description="Count of variant callers which listed this variant in the Jacquard tagged VCF">'.format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_REPORTED),
                          split_metaheader[0])

    def test_add_tag_values(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:1:1|Y:1:1\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._CallersReportedTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}{}|X:1:1:2|Y:1:1:2\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED, summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_REPORTED))
        self.assertEquals(expected, processedVcfRecord.text())

    def test_add_tag_values_nullValues(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:.:.|Y:.:.\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._CallersReportedTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}{}|X:.:.:0|Y:.:.:0\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED, summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_REPORTED))
        self.assertEquals(expected, processedVcfRecord.text())

class CallersPassedListTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = summarize_caller._CallersPassedListTag().metaheader.split("\n")
        self.assertEquals('##FORMAT=<ID={}CALLERS_PASSED_LIST,Number=.,Type=String,Description="Comma-separated list of variant caller short-names where FILTER = PASS for this variant in the Jacquard tagged VCF">'.format(summarize_caller.JQ_SUMMARY_TAG),
                          split_metaheader[0])

    def test_add_tag_values(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:1:1|Y:1:0\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._CallersPassedListTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}CALLERS_PASSED_LIST|X:1:1:MT,VS|Y:1:0:MT\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED, summarize_caller.JQ_SUMMARY_TAG))
        self.assertEquals(expected, processedVcfRecord.text())

    def test_add_tag_values_NoCallersPassed(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:0:0|Y:0:0\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._CallersPassedListTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}CALLERS_PASSED_LIST|X:0:0:.|Y:0:0:.\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED, summarize_caller.JQ_SUMMARY_TAG))
        self.assertEquals(expected, processedVcfRecord.text())

    def test_add_tag_values_nullValues(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:1:.|Y:1:.\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._CallersPassedListTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}CALLERS_PASSED_LIST|X:1:.:MT|Y:1:.:MT\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED, summarize_caller.JQ_SUMMARY_TAG))
        self.assertEquals(expected, processedVcfRecord.text())

class CallersPassedTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = summarize_caller._CallersPassedTag().metaheader.split("\n")
        self.assertEquals('##FORMAT=<ID={}{},Number=1,Type=Integer,Description="Count of variant callers where FILTER = PASS for this variant in the Jacquard tagged VCF">'.format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_PASSED),
                          split_metaheader[0])

    def test_add_tag_values(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:1:1|Y:1:0\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._CallersPassedTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}{}|X:1:1:2|Y:1:0:1\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED, summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_PASSED))
        self.assertEquals(expected, processedVcfRecord.text())

    def test_add_tag_values_nullValues(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:.:.|Y:.:.\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._CallersPassedTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}{}|X:.:.:0|Y:.:.:0\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED, summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_PASSED))
        self.assertEquals(expected, processedVcfRecord.text())

class SamplesReportedTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = summarize_caller._SamplesReported().metaheader.split("\n")
        self.assertEquals('##INFO=<ID={}{},Number=1,Type=Integer,Description="Count of samples where this variant appeared in any of the Jacquard tagged VCFs (regardless of quality/filtering)">'.format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_SAMPLES_REPORTED),
                          split_metaheader[0])

    def test_add_tag_values(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}|X:2|Y:1\n".format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_REPORTED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._SamplesReported()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO;{}{}=2|JQ_DP:{}{}|X:2|Y:1\n".format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_SAMPLES_REPORTED, summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_REPORTED))
        self.assertEquals(expected, processedVcfRecord.text())

    def test_add_tag_values_nullValues(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}|X:.|Y:.\n".format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_REPORTED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._SamplesReported()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO;{}{}=0|JQ_DP:{}{}|X:.|Y:.\n".format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_SAMPLES_REPORTED, summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_REPORTED))
        self.assertEquals(expected, processedVcfRecord.text())

class SamplesPassedTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = summarize_caller._SamplesPassed().metaheader.split("\n")
        self.assertEquals('##INFO=<ID={}{},Number=1,Type=Integer,Description="Count of samples where a variant caller passed the filter in any of the Jacquard tagged VCFs">'.format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_SAMPLES_PASSED),
                          split_metaheader[0])

    def test_add_tag_values_onePassed(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}|X:2|Y:0\n".format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_PASSED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._SamplesPassed()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO;{}{}=1|JQ_DP:{}{}|X:2|Y:0\n".format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_SAMPLES_PASSED, summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_PASSED))
        self.assertEquals(expected, processedVcfRecord.text())

    def test_add_tag_values_nonePassed(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}|X:0|Y:0\n".format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_PASSED))
        vcf_record = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._SamplesPassed()
        tag.add_tag_values(vcf_record)

        info_tag = summarize_caller.JQ_SUMMARY_TAG + summarize_caller.JQ_SAMPLES_PASSED
        self.assertIn(info_tag, vcf_record.info_dict)
        self.assertEquals("0", vcf_record.info_dict[info_tag])

    def test_add_tag_values_nullValues(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}|X:.|Y:.\n".format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_PASSED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = summarize_caller._SamplesPassed()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO;{}{}=0|JQ_DP:{}{}|X:.|Y:.\n".format(summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_SAMPLES_PASSED, summarize_caller.JQ_SUMMARY_TAG, summarize_caller.JQ_PASSED))
        self.assertEquals(expected, processedVcfRecord.text())

class HCGenotypeTagTestCase(test_case.JacquardBaseTestCase):
    _TAG_ID = "{}HC_GT".format(summarize_caller.JQ_SUMMARY_TAG)

    def test_metaheader(self):
        self.assertEquals('##FORMAT=<ID={}HC_GT,Number=1,Type=String,Description="High confidence consensus genotype (inferred from JQ_*_GT and JQ_*_CALLER_PASSED). Majority rules; ties go to the least unusual variant (0/1>0/2>1/1). Variants which failed their filter are ignored.">'.format(summarize_caller.JQ_SUMMARY_TAG),
                          summarize_caller._HCGenotypeTag().metaheader)

    def test_add_tag_values(self):
        tag = summarize_caller._HCGenotypeTag()
        sample_tag_values = {"SA": {"JQ_foo_GT": "0/1", "JQ_bar_GT": "0/1", "JQ_baz_GT":"."},
                             "SB": {"JQ_foo_GT": "0/0", "JQ_bar_GT": "0/0", "JQ_baz_GT":"0/1"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals("0/1", record.sample_tag_values["SA"][HCGenotypeTagTestCase._TAG_ID])
        self.assertEquals("0/0", record.sample_tag_values["SB"][HCGenotypeTagTestCase._TAG_ID])

    def test_add_tag_values_allNulls(self):
        tag = summarize_caller._HCGenotypeTag()
        sample_tag_values = {"SA": {"JQ_foo_GT": ".", "JQ_bar_GT": ".", "JQ_baz_GT":"."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals(".", record.sample_tag_values["SA"][HCGenotypeTagTestCase._TAG_ID])

    def test_add_tag_values_onlyCountGTTags(self):
        tag = summarize_caller._HCGenotypeTag()
        sample_tag_values = {"SA": {"JQ_foo_XX": "0/1", "JQ_bar_GT": ".", "JQ_baz_GT":"0/0"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("0/0", record.sample_tag_values["SA"][HCGenotypeTagTestCase._TAG_ID])

    def test_add_tag_values_missingGTTags(self):
        tag = summarize_caller._HCGenotypeTag()
        sample_tag_values = {"SA": {"JQ_foo_XX": "0/1"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals(".", record.sample_tag_values["SA"][HCGenotypeTagTestCase._TAG_ID])

    def test_prioritize_genotype(self):
        tag = summarize_caller._HCGenotypeTag()
        self.assertEquals("0/1", tag._prioritize_genotype(["0/1"]))
        self.assertEquals("0/1", tag._prioritize_genotype(["0/0", "0/1"]))
        self.assertEquals("1/1", tag._prioritize_genotype(["0/0", "1/1"]))
        self.assertEquals("1/1", tag._prioritize_genotype(["1/1", "1/1", "0/1"]))
        self.assertEquals("0/1", tag._prioritize_genotype(["1/1", "0/1", "0/1"]))
        self.assertEquals("0/1", tag._prioritize_genotype(["0/1", "0/2"]))

class AlleleFreqRangeTagTestCase(test_case.JacquardBaseTestCase):
    _TAG_ID = "{}AF_RANGE".format(summarize_caller.JQ_SUMMARY_TAG)

    def test_metaheader(self):
        split_meta_header = summarize_caller._AlleleFreqRangeTag().metaheader.split("\n")
        self.assertEqual('##FORMAT=<ID={0}AF_RANGE,Number=1,Type=Float,Description="Max(allele frequency) - min (allele frequency) across recognized callers.">'.format(summarize_caller.JQ_SUMMARY_TAG),
                         split_meta_header[0])

    def test_add_tag_values(self):
        tag = summarize_caller._AlleleFreqRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": "0", "JQ_bar_AF": "0.1", "JQ_baz_AF":"0.2"},
                             "SB": {"JQ_foo_AF": "0.2", "JQ_bar_AF": "0.3", "JQ_baz_AF":"0.4"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("0.2", record.sample_tag_values["SA"][AlleleFreqRangeTagTestCase._TAG_ID])
        self.assertEquals("0.2", record.sample_tag_values["SB"][AlleleFreqRangeTagTestCase._TAG_ID])

    def test_add_tag_values_nullsDoNotCount(self):
        tag = summarize_caller._AlleleFreqRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": ".", "JQ_bar_AF": "0.1", "JQ_baz_AF":"0.2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals("0.1", record.sample_tag_values["SA"][AlleleFreqRangeTagTestCase._TAG_ID])

    def test_add_tag_values_oneValueReturnsNull(self):
        tag = summarize_caller._AlleleFreqRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": "."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals(".", record.sample_tag_values["SA"][AlleleFreqRangeTagTestCase._TAG_ID])

    def test_add_tag_values_OnlyCountAFTags(self):
        tag = summarize_caller._AlleleFreqRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_XX": "0", "JQ_bar_AF": "0.1", "JQ_baz_AF":"0.2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("0.1", record.sample_tag_values["SA"][AlleleFreqRangeTagTestCase._TAG_ID])

    def test_add_tag_values_allNulls(self):
        tag = summarize_caller._AlleleFreqRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": ".", "JQ_bar_AF": ".", "JQ_baz_AF":"."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals(".", record.sample_tag_values["SA"][AlleleFreqRangeTagTestCase._TAG_ID])

    def test_add_tag_values_missingTag(self):
        tag = summarize_caller._AlleleFreqRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_XX": "."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals(".", record.sample_tag_values["SA"][AlleleFreqRangeTagTestCase._TAG_ID])

    def test_add_tag_values_multAlts(self):
        tag = summarize_caller._AlleleFreqRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": "0,0.1", "JQ_bar_AF": "0.1,0.2"},
                             "SB": {"JQ_foo_AF": "0.2", "JQ_bar_AF": "0.3"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("0.1,0.1", record.sample_tag_values["SA"][AlleleFreqRangeTagTestCase._TAG_ID])

    def test_add_tag_values_rounds(self):
        tag = summarize_caller._AlleleFreqRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": "0", "JQ_bar_AF": "0.666666"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("0.67", record.sample_tag_values["SA"][AlleleFreqRangeTagTestCase._TAG_ID])

    def test_add_tag_values_inconsistentMultAlt(self):
        tag = summarize_caller._AlleleFreqRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": "0,0.1", "JQ_bar_AF": "0.2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        self.assertRaisesRegexp(JQException,
                                r"Error summarizing values \[.*\] at record \[.*\]",
                                tag.add_tag_values,
                                record)

class AlleleFreqAverageTagTestCase(test_case.JacquardBaseTestCase):
    _TAG_ID = "{}AF_AVERAGE".format(summarize_caller.JQ_SUMMARY_TAG)

    def test_metaheader(self):
        split_meta_header = summarize_caller._AlleleFreqAverageTag().metaheader.split("\n")
        self.assertEqual('##FORMAT=<ID={0}AF_AVERAGE,Number=1,Type=Float,Description="Average allele frequency across recognized variant callers that reported frequency for this position [average(JQ_*_AF)].">'.format(summarize_caller.JQ_SUMMARY_TAG),
                         split_meta_header[0])

    def test_add_tag_values(self):
        tag = summarize_caller._AlleleFreqAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": "0", "JQ_bar_AF": "0.1", "JQ_baz_AF":"0.2"},
                             "SB": {"JQ_foo_AF": "0.2", "JQ_bar_AF": "0.3", "JQ_baz_AF":"0.4"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("0.1", record.sample_tag_values["SA"][AlleleFreqAverageTagTestCase._TAG_ID])
        self.assertEquals("0.3", record.sample_tag_values["SB"][AlleleFreqAverageTagTestCase._TAG_ID])

    def test_add_tag_values_nullsDoNotCount(self):
        tag = summarize_caller._AlleleFreqAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": ".", "JQ_bar_AF": ".", "JQ_baz_AF":"0.2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("0.2", record.sample_tag_values["SA"][AlleleFreqAverageTagTestCase._TAG_ID])

    def test_add_tag_values_OnlyCountAFTags(self):
        tag = summarize_caller._AlleleFreqAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_XX": "0", "JQ_bar_XX": "0.1", "JQ_baz_AF":"0.2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("0.2", record.sample_tag_values["SA"][AlleleFreqAverageTagTestCase._TAG_ID])

    def test_add_tag_values_allNulls(self):
        tag = summarize_caller._AlleleFreqAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": ".", "JQ_bar_AF": ".", "JQ_baz_AF":"."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals(".", record.sample_tag_values["SA"][AlleleFreqAverageTagTestCase._TAG_ID])

    def test_add_tag_values_missingTag(self):
        tag = summarize_caller._AlleleFreqAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_XX": "."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals(".", record.sample_tag_values["SA"][AlleleFreqAverageTagTestCase._TAG_ID])

    def test_add_tag_values_multAlts(self):
        tag = summarize_caller._AlleleFreqAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": "0,0.1", "JQ_bar_AF": "0.1,0.2"},
                             "SB": {"JQ_foo_AF": "0.2", "JQ_bar_AF": "0.3"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("0.05,0.15", record.sample_tag_values["SA"][AlleleFreqAverageTagTestCase._TAG_ID])

    def test_add_tag_values_rounds(self):
        tag = summarize_caller._AlleleFreqAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": "0", "JQ_bar_AF": "0.666666"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("0.33", record.sample_tag_values["SA"][AlleleFreqAverageTagTestCase._TAG_ID])

    def test_add_tag_values_inconsistentMultAlt(self):
        tag = summarize_caller._AlleleFreqAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_AF": "0,0.1", "JQ_bar_AF": "0.2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        self.assertRaisesRegexp(JQException,
                                r"Error summarizing values \[.*\] at record \[.*\]",
                                tag.add_tag_values,
                                record)

class DepthRangeTagTestCase(test_case.JacquardBaseTestCase):
    _TAG_ID = "{}DP_RANGE".format(summarize_caller.JQ_SUMMARY_TAG)

    def test_metaheader(self):
        split_meta_header = summarize_caller._DepthRangeTag().metaheader.split("\n")
        self.assertEquals('##FORMAT=<ID={0}DP_RANGE,Number=1,Type=Float,Description="Max(depth) - min (depth) across recognized callers.">'.format(summarize_caller.JQ_SUMMARY_TAG),
                          split_meta_header[0])

    def test_add_tag_values(self):
        tag = summarize_caller._DepthRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": "0", "JQ_bar_DP": "1", "JQ_baz_DP":"2"},
                             "SB": {"JQ_foo_DP": "2", "JQ_bar_DP": "3", "JQ_baz_DP":"4"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("2", record.sample_tag_values["SA"][DepthRangeTagTestCase._TAG_ID])
        self.assertEquals("2", record.sample_tag_values["SB"][DepthRangeTagTestCase._TAG_ID])

    def test_add_tag_values_oneValueReturnsNull(self):
        tag = summarize_caller._DepthRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": "."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals(".", record.sample_tag_values["SA"][DepthRangeTagTestCase._TAG_ID])

    def test_add_tag_values_nullsDoNotCount(self):
        tag = summarize_caller._DepthRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": ".", "JQ_bar_DP": "1", "JQ_baz_DP":"2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("1", record.sample_tag_values["SA"][DepthRangeTagTestCase._TAG_ID])

    def test_add_tag_values_OnlyCountDPTags(self):
        tag = summarize_caller._DepthRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_XX": "0", "JQ_bar_DP": "1", "JQ_baz_DP":"2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("1", record.sample_tag_values["SA"][DepthRangeTagTestCase._TAG_ID])

    def test_add_tag_values_allNulls(self):
        tag = summarize_caller._DepthRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": ".", "JQ_bar_DP": ".", "JQ_baz_DP":"."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals(".", record.sample_tag_values["SA"][DepthRangeTagTestCase._TAG_ID])

    def test_add_tag_values_missingTag(self):
        tag = summarize_caller._DepthRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_XX": "."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals(".", record.sample_tag_values["SA"][DepthRangeTagTestCase._TAG_ID])

    def test_add_tag_values_multAlts(self):
        tag = summarize_caller._DepthRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": "0,1", "JQ_bar_DP": "2,3"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("2,2", record.sample_tag_values["SA"][DepthRangeTagTestCase._TAG_ID])

    def test_add_tag_values_rounds(self):
        tag = summarize_caller._DepthRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": "0", "JQ_bar_DP": "42.666666"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("42.67", record.sample_tag_values["SA"][DepthRangeTagTestCase._TAG_ID])

    def test_add_tag_values_inconsistentMultAlt(self):
        tag = summarize_caller._DepthRangeTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": "0,0.1", "JQ_bar_DP": "0.2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        self.assertRaisesRegexp(JQException,
                                r"Error summarizing values \[.*\] at record \[.*\]",
                                tag.add_tag_values,
                                record)

class DepthAverageTagTestCase(test_case.JacquardBaseTestCase):
    _TAG_ID = "{}DP_AVERAGE".format(summarize_caller.JQ_SUMMARY_TAG)

    def test_metaheader(self):
        split_meta_header = summarize_caller._DepthAverageTag().metaheader.split("\n")
        self.assertEquals('##FORMAT=<ID={}DP_AVERAGE,Number=1,Type=Float,Description="Average allele frequency across recognized variant callers that reported frequency for this position; rounded to integer [round(average(JQ_*_DP))].">'.format(summarize_caller.JQ_SUMMARY_TAG),
                          split_meta_header[0])

    def test_add_tag_values(self):
        tag = summarize_caller._DepthAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": "0", "JQ_bar_DP": "1", "JQ_baz_DP":"2"},
                             "SB": {"JQ_foo_DP": "2", "JQ_bar_DP": "3", "JQ_baz_DP":"4"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("1", record.sample_tag_values["SA"][DepthAverageTagTestCase._TAG_ID])
        self.assertEquals("3", record.sample_tag_values["SB"][DepthAverageTagTestCase._TAG_ID])

    def test_add_tag_values_oneValueReturnsNull(self):
        tag = summarize_caller._DepthAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": "."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals(".", record.sample_tag_values["SA"][DepthAverageTagTestCase._TAG_ID])

    def test_add_tag_values_nullsDoNotCount(self):
        tag = summarize_caller._DepthAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": ".", "JQ_bar_DP": "1", "JQ_baz_DP":"2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("1.5", record.sample_tag_values["SA"][DepthAverageTagTestCase._TAG_ID])

    def test_add_tag_values_OnlyCountDPTags(self):
        tag = summarize_caller._DepthAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_XX": "0", "JQ_bar_DP": "1", "JQ_baz_DP":"2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("1.5", record.sample_tag_values["SA"][DepthAverageTagTestCase._TAG_ID])

    def test_add_tag_values_allNulls(self):
        tag = summarize_caller._DepthAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": ".", "JQ_bar_DP": ".", "JQ_baz_DP":"."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals(".", record.sample_tag_values["SA"][DepthAverageTagTestCase._TAG_ID])

    def test_add_tag_values_missingTag(self):
        tag = summarize_caller._DepthAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_XX": "."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals(".", record.sample_tag_values["SA"][DepthAverageTagTestCase._TAG_ID])

    def test_add_tag_values_multAlts(self):
        tag = summarize_caller._DepthAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": "0,1", "JQ_bar_DP": "2,3"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("1,2", record.sample_tag_values["SA"][DepthAverageTagTestCase._TAG_ID])

    def test_add_tag_values_rounds(self):
        tag = summarize_caller._DepthAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": "0", "JQ_bar_DP": "3"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("1.5", record.sample_tag_values["SA"][DepthAverageTagTestCase._TAG_ID])

    def test_add_tag_values_inconsistentMultAlt(self):
        tag = summarize_caller._DepthAverageTag()
        sample_tag_values = {"SA": {"JQ_foo_DP": "0,0.1", "JQ_bar_DP": "0.2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        self.assertRaisesRegexp(JQException,
                                r"Error summarizing values \[.*\] at record \[.*\]",
                                tag.add_tag_values,
                                record)

class SomaticTagTestCase(test_case.JacquardBaseTestCase):
    _TAG_ID = "{}SOM_COUNT".format(summarize_caller.JQ_SUMMARY_TAG)

    def test_metaheader(self):
        split_meta_header = summarize_caller._SomaticTag().metaheader.split("\n")
        self.assertEqual('##FORMAT=<ID={0}SOM_COUNT,Number=1,Type=Integer,' \
                      'Description="Count of recognized variant callers ' \
                      'that reported confident somatic call for this '\
                      'sample-position.">'\
                      .format(summarize_caller.JQ_SUMMARY_TAG), split_meta_header[0])

    def test_add_tag_values(self):
        tag = summarize_caller._SomaticTag()
        sample_tag_values = {"SA": {"JQ_foo_HC_SOM": "1", "JQ_bar_HC_SOM": "1"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals("2", record.sample_tag_values["SA"][SomaticTagTestCase._TAG_ID])

    def test_add_tag_values_allZero(self):
        tag = summarize_caller._SomaticTag()
        sample_tag_values = {"SA": {"JQ_foo_HC_SOM": "0", "JQ_bar_HC_SOM": "0"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals("0", record.sample_tag_values["SA"][SomaticTagTestCase._TAG_ID])

    def test_add_tag_values_nullsDoNotCount(self):
        tag = summarize_caller._SomaticTag()
        sample_tag_values = {"SA": {"JQ_foo_HC_SOM": "1", "JQ_bar_HC_SOM": "."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals("1", record.sample_tag_values["SA"][SomaticTagTestCase._TAG_ID])

    def test_add_tag_values_OnlyCountHCSOMTags(self):
        tag = summarize_caller._SomaticTag()
        sample_tag_values = {"SA": {"JQ_foo": "1", "JQ_foo_HC_SOM": "1", "JQ_bar_HC_SOM": "1"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals("2", record.sample_tag_values["SA"][SomaticTagTestCase._TAG_ID])

    def test_add_tag_values_allNulls(self):
        tag = summarize_caller._SomaticTag()
        sample_tag_values = {"SA": {"JQ_foo_HC_SOM": ".", "JQ_bar_HC_SOM": "."}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals(".", record.sample_tag_values["SA"][SomaticTagTestCase._TAG_ID])

    def test_add_tag_values_missingTag(self):
        tag = summarize_caller._SomaticTag()
        sample_tag_values = {"SA": {"JQ_foo": "1",}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)
        self.assertEquals(".", record.sample_tag_values["SA"][SomaticTagTestCase._TAG_ID])

    def test_add_tag_values_multAlts(self):
        tag = summarize_caller._SomaticTag()
        sample_tag_values = {"SA": {"JQ_foo_HC_SOM": "0,1", "JQ_bar_HC_SOM": "1,1"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        tag.add_tag_values(record)

        self.assertEquals("1,2", record.sample_tag_values["SA"][SomaticTagTestCase._TAG_ID])

    def test_add_tag_values_inconsistentMultAlt(self):
        tag = summarize_caller._SomaticTag()
        sample_tag_values = {"SA": {"JQ_foo_HC_SOM": "0,0.1", "JQ_bar_HC_SOM": "0.2"}}
        record = VcfRecord("CHROM", "POS", "REF", "ALT",
                           sample_tag_values=sample_tag_values)
        self.assertRaisesRegexp(JQException,
                                r"Error summarizing values \[.*\] at record \[.*\]",
                                tag.add_tag_values,
                                record)

class SummarizeCallerTestCase(test_case.JacquardBaseTestCase):
    @staticmethod
    def _sum_to_string(collection_of_numbers):
        return str(sum(collection_of_numbers))

    def test_aggregate_numeric_values_ints(self):
        input_values = ["0", "1", "2", "4"]
        actual_value = summarize_caller._aggregate_numeric_values(input_values,
                                                                  self._sum_to_string)
        self.assertEquals("7", actual_value)

    def test_aggregate_numeric_values_floats(self):
        input_values = ["0.0", "0.1", "0.2", "0.4"]
        actual_value = summarize_caller._aggregate_numeric_values(input_values,
                                                                  self._sum_to_string)
        self.assertEquals("0.7", actual_value)

    def test_aggregate_numeric_values_listsOfInts(self):
        input_values = ["0,1", "2,4"]
        actual_value = summarize_caller._aggregate_numeric_values(input_values,
                                                                  self._sum_to_string)
        self.assertEquals("2,5", actual_value)

    def test_aggregate_numeric_values_listsOfFloats(self):
        input_values = ["0.0,0.1", "0.2,0.4"]
        actual_value = summarize_caller._aggregate_numeric_values(input_values,
                                                                  self._sum_to_string)
        self.assertEquals("0.2,0.5", actual_value)

    def test_aggregate_numeric_values_listsOf3(self):
        input_values = ["0,1,3", "4,5,6"]
        actual_value = summarize_caller._aggregate_numeric_values(input_values,
                                                                  self._sum_to_string)
        self.assertEquals("4,6,9", actual_value)

    def test_get_non_null_values(self):
        sample_tag_values = {"SA": {"JQ_A_AF":"0",
                                    "JQ_B_AF":"1",
                                    "JQ_C_AF":"2"}}
        record = VcfRecord("chr1", "42", "A", "C", sample_tag_values=sample_tag_values)
        actual_values = summarize_caller._get_non_null_values(record,
                                                              "SA",
                                                              common_tags.ALLELE_FREQ_TAG)
        self.assertEquals(["0", "1", "2"], sorted(actual_values))

    def test_get_non_null_values_hasNulls(self):
        sample_tag_values = {"SA": {"JQ_A_AF":".",
                                    "JQ_B_AF":".",
                                    "JQ_C_AF":"2"}}
        record = VcfRecord("chr1", "42", "A", "C", sample_tag_values=sample_tag_values)
        actual_values = summarize_caller._get_non_null_values(record,
                                                              "SA",
                                                              common_tags.ALLELE_FREQ_TAG)
        self.assertEquals(["2"], actual_values)

    def test_get_non_null_values_invalidSample(self):
        sample_tag_values = {"SA": {"JQ_A_AF":".",
                                    "JQ_B_AF":".",
                                    "JQ_C_AF":"2"}}
        record = VcfRecord("chr1", "42", "A", "C", sample_tag_values=sample_tag_values)
        self.assertRaisesRegexp(utils.JQException,
                                r"Sample \[FOO\] was not recognized",
                                summarize_caller._get_non_null_values,
                                record,
                                "FOO",
                                re.compile("^JQ_.*_AF$"))

    def test_get_non_null_values_missingTag(self):
        sample_tag_values = {"SA": {"JQ_A_AF":".",
                                    "JQ_B_AF":".",
                                    "JQ_C_AF":"2"}}
        record = VcfRecord("chr1", "42", "A", "C", sample_tag_values=sample_tag_values)
        actual_values = summarize_caller._get_non_null_values(record,
                                                              "SA",
                                                              common_tags.DEPTH_TAG)
        self.assertEquals([], actual_values)

    def test_summary_tag_prefix(self):
        self.assertEquals("JQ_SUMMARY_", summarize_caller.JQ_SUMMARY_TAG)

    def test_add_tags(self):
        sample_tag_values = {"SA": {"JQ_foo_AF":"0", "JQ_VS_CALLER_REPORTED":"1", 'JQ_MT_CALLER_REPORTED':"1", 'JQ_VS_CALLER_PASSED': "0", 'JQ_MT_CALLER_PASSED': "0"},
                             "SB": {"JQ_foo_AF":"0.2", "JQ_VS_CALLER_REPORTED":"1", 'JQ_MT_CALLER_REPORTED':"1", 'JQ_VS_CALLER_PASSED': "0", 'JQ_MT_CALLER_PASSED': "0"}}
        record = VcfRecord("chr1", "42", "A", "C", sample_tag_values=sample_tag_values)
        actual_record = summarize_caller.SummarizeCaller().add_tags(record)

#         self.assertEquals("2", actual_record.info_dict["JQ_SUMMARY_SAMPLES_REPORTED_COUNT"])
#         self.assertEquals("0", actual_record.info_dict["JQ_SUMMARY_SAMPLES_PASSED_COUNT"])

        self.assertEquals("2", actual_record.sample_tag_values["SA"]["JQ_SUMMARY_CALLERS_REPORTED_COUNT"])
        self.assertEquals("0", actual_record.sample_tag_values["SA"]["JQ_SUMMARY_CALLERS_PASSED_COUNT"])
        self.assertEquals("0", actual_record.sample_tag_values["SA"]["JQ_SUMMARY_AF_AVERAGE"])
        self.assertEquals(".", actual_record.sample_tag_values["SA"]["JQ_SUMMARY_AF_RANGE"])
        self.assertEquals(".", actual_record.sample_tag_values["SA"]["JQ_SUMMARY_DP_AVERAGE"])
        self.assertEquals(".", actual_record.sample_tag_values["SA"]["JQ_SUMMARY_DP_RANGE"])
        self.assertEquals(".", actual_record.sample_tag_values["SA"]["JQ_SUMMARY_SOM_COUNT"])

        self.assertEquals("2", actual_record.sample_tag_values["SB"]["JQ_SUMMARY_CALLERS_REPORTED_COUNT"])
        self.assertEquals("0", actual_record.sample_tag_values["SB"]["JQ_SUMMARY_CALLERS_PASSED_COUNT"])
        self.assertEquals("0.2", actual_record.sample_tag_values["SB"]["JQ_SUMMARY_AF_AVERAGE"])
        self.assertEquals(".", actual_record.sample_tag_values["SB"]["JQ_SUMMARY_AF_RANGE"])
        self.assertEquals(".", actual_record.sample_tag_values["SB"]["JQ_SUMMARY_DP_AVERAGE"])
        self.assertEquals(".", actual_record.sample_tag_values["SB"]["JQ_SUMMARY_DP_RANGE"])
        self.assertEquals(".", actual_record.sample_tag_values["SB"]["JQ_SUMMARY_SOM_COUNT"])

    def test_get_new_metaheaders(self):
        expected = ('##FORMAT=<ID={}{},'
                    'Number=1,'
                    'Type=Integer,'
                    'Description="Count of variant callers which listed this variant in the Jacquard tagged VCF">')\
                    .format(summarize_caller.JQ_SUMMARY_TAG,
                            summarize_caller.JQ_REPORTED)

        actual = summarize_caller.SummarizeCaller().get_metaheaders()

        split_actual = actual[0].split("\n")
        first_meta_header = split_actual[0]

        self.assertEqual(expected, first_meta_header)
        self.assertEqual(12, len(actual))
        self.assertEqual(1, len(split_actual))


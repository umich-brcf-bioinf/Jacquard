#pylint: disable=too-few-public-methods, invalid-name, line-too-long
#pylint: disable=too-many-instance-attributes, too-many-public-methods
from __future__ import print_function, absolute_import

from collections import OrderedDict
import os
import unittest

import jacquard.variant_callers.common_tags as common_tags
import jacquard.variant_callers.jacquard_consensus_caller as consensus_caller
import jacquard.variant_callers.mutect as mutect
import jacquard.variant_callers.varscan as varscan
from jacquard.vcf import VcfRecord
from jacquard import __version__
import test.test_case as test_case


class MockVcfRecord(object):
    def __init__(self, content):
        content = content.split("\t")
        self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, \
            self.filter, self.info, self.format = content[0:9]
        self.samples = content[9:]

        tags = self.format.split(":")
        self.format_set = tags

        self.sample_dict = {}
        for i, sample in enumerate(self.samples):
            values = sample.split(":")
            self.sample_dict[i] = OrderedDict(zip(tags, values))

    def insert_format_field(self, fieldname, field_dict):
        if fieldname in self.format_set:
            raise KeyError
        self.format_set.append(fieldname)

        if field_dict.keys() != self.sample_dict.keys():
            raise KeyError()
        for key in self.sample_dict.keys():
            self.sample_dict[key][fieldname] = str(field_dict[key])

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

class CallersReportedListTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = consensus_caller._CallersReportedListTag().metaheader.split("\n")
        self.assertEquals('##FORMAT=<ID={}{},Number=.,Type=String,Description="Comma-separated list variant callers which listed this variant in the Jacquard tagged VCF",Source="Jacquard",Version="{}">'.format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_REPORTED_LIST, __version__),
                          split_metaheader[0])

    def test_add_tag_values(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:1:1|Y:1:1\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED_TAG, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED_TAG))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = consensus_caller._CallersReportedListTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}{}|X:1:1:MT,VS|Y:1:1:MT,VS\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED_TAG, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED_TAG, consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_REPORTED_LIST))
        self.assertEquals(expected, processedVcfRecord.asText())

class CallersReportedTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = consensus_caller._CallersReportedTag().metaheader.split("\n")
        self.assertEquals('##FORMAT=<ID={}{},Number=1,Type=Integer,Description="Count of variant callers which listed this variant in the Jacquard tagged VCF",Source="Jacquard",Version="{}">'.format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_REPORTED, __version__),
                          split_metaheader[0])

    def test_add_tag_values(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:1:1|Y:1:1\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED_TAG, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED_TAG))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = consensus_caller._CallersReportedTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}{}|X:1:1:2|Y:1:1:2\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_REPORTED_TAG, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_REPORTED_TAG, consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_REPORTED))
        self.assertEquals(expected, processedVcfRecord.asText())

class CallersPassedListTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = consensus_caller._CallersPassedListTag().metaheader.split("\n")
        self.assertEquals('##FORMAT=<ID={}{},Number=.,Type=String,Description="Comma-separated list of variant caller short-names where FILTER = PASS for this variant in the Jacquard tagged VCF",Source="Jacquard",Version="{}">'.format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_PASSED_LIST, __version__),
                          split_metaheader[0])

    def test_add_tag_values(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:1:1|Y:1:0\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED_TAG, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED_TAG))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = consensus_caller._CallersPassedListTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}{}|X:1:1:MT,VS|Y:1:0:MT\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED_TAG, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED_TAG, consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_PASSED_LIST))
        self.assertEquals(expected, processedVcfRecord.asText())

class CallersPassedTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = consensus_caller._CallersPassedTag().metaheader.split("\n")
        self.assertEquals('##FORMAT=<ID={}{},Number=1,Type=Integer,Description="Count of variant callers where FILTER = PASS for this variant in the Jacquard tagged VCF",Source="Jacquard",Version="{}">'.format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_PASSED, __version__),
                          split_metaheader[0])

    def test_add_tag_values(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}|X:1:1|Y:1:0\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED_TAG, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED_TAG))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = consensus_caller._CallersPassedTag()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}:{}{}:{}{}|X:1:1:2|Y:1:0:1\n".format(mutect.JQ_MUTECT_TAG, common_tags.CALLER_PASSED_TAG, varscan.JQ_VARSCAN_TAG, common_tags.CALLER_PASSED_TAG, consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_PASSED))
        self.assertEquals(expected, processedVcfRecord.asText())

class SamplesReportedTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = consensus_caller._SamplesReported().metaheader.split("\n")
        self.assertEquals('##INFO=<ID={}{},Number=1,Type=Integer,Description="Count of samples where this variant appeared in any of the Jacquard tagged VCFs (regardless of quality/filtering)",Source="Jacquard",Version="{}">'.format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_SAMPLES_REPORTED, __version__),
                          split_metaheader[0])

    def test_add_tag_values(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}|X:2|Y:1\n".format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_REPORTED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = consensus_caller._SamplesReported()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO;{}{}=2|JQ_DP:{}{}|X:2|Y:1\n".format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_SAMPLES_REPORTED, consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_REPORTED))
        self.assertEquals(expected, processedVcfRecord.asText())

class SamplesPassedTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_metaheader = consensus_caller._SamplesPassed().metaheader.split("\n")
        self.assertEquals('##INFO=<ID={}{},Number=1,Type=Integer,Description="Count of samples where a variant caller passed the filter in any of the Jacquard tagged VCFs",Source="Jacquard",Version="{}">'.format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_SAMPLES_PASSED, __version__),
                          split_metaheader[0])

    def test_add_tag_values_onePassed(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}|X:2|Y:0\n".format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_PASSED))
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = consensus_caller._SamplesPassed()
        tag.add_tag_values(processedVcfRecord)

        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO;{}{}=1|JQ_DP:{}{}|X:2|Y:0\n".format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_SAMPLES_PASSED, consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_PASSED))
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_add_tag_values_nonePassed(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{}{}|X:0|Y:0\n".format(consensus_caller.JQ_CONSENSUS_TAG, consensus_caller.JQ_PASSED))
        vcf_record = VcfRecord.parse_record(line, ["SA", "SB"])
        tag = consensus_caller._SamplesPassed()
        tag.add_tag_values(vcf_record)

        info_tag = consensus_caller.JQ_CONSENSUS_TAG + consensus_caller.JQ_SAMPLES_PASSED
        self.assertIn(info_tag, vcf_record.info_dict)
        self.assertEquals("0", vcf_record.info_dict[info_tag])


class AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        split_meta_header = consensus_caller._AlleleFreqTag().metaheader.split("\n")
        self.assertEqual('##FORMAT=<ID={0}AF_AVERAGE,Number=1,Type=Float,Description="Average allele frequency across recognized variant callers that reported frequency for this position [average(JQ_*_AF)].",Source="Jacquard",Version="{1}">'.format(consensus_caller.JQ_CONSENSUS_TAG, __version__),
                         split_meta_header[0])

    def test_insert_consensus(self):
        tag = consensus_caller._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF|X:0:0.1:0.2|Y:0.2:0.3:0.4\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF:{0}AF_AVERAGE:{0}AF_RANGE|X:0:0.1:0.2:0.1:0.2|Y:0.2:0.3:0.4:0.3:0.2\n".format(consensus_caller.JQ_CONSENSUS_TAG).replace('|', "\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_multAlts(self):
        tag = consensus_caller._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF|0,0:0.2,0.4|0,0:0.1,0.3\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF:{0}AF_AVERAGE:{0}AF_RANGE|0,0:0.2,0.4:0.1,0.2:0.2,0.4|0,0:0.1,0.3:0.05,0.15:0.1,0.3\n".format(consensus_caller.JQ_CONSENSUS_TAG).replace('|', "\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_noJQ_AFTags(self):
        tag = consensus_caller._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP|X|Y\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{0}AF_AVERAGE:{0}AF_RANGE|X:.:.|Y:.:.\n".format(consensus_caller.JQ_CONSENSUS_TAG).replace('|', "\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_noNullValuesInAvgCalculation(self):
        tag = consensus_caller._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF|X:0:0.1:.|Y:0.2:0.3:.\n".replace('|', "\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF:{0}AF_AVERAGE:{0}AF_RANGE|X:0:0.1:.:0.05:0.1|Y:0.2:0.3:.:0.25:0.1\n".format(consensus_caller.JQ_CONSENSUS_TAG).replace('|', "\t")
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class DepthTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_meta_header = consensus_caller._DepthTag().metaheader.split("\n")
        self.assertEqual('##FORMAT=<ID={0}DP_AVERAGE,Number=1,Type=Float,' \
                      'Description="Average allele frequency across ' \
                      'recognized variant callers that reported ' \
                      'frequency for this position; rounded to integer '\
                      '[round(average(JQ_*_DP))].",' \
                      'Source="Jacquard",Version="{1}">'\
                      .format(consensus_caller.JQ_CONSENSUS_TAG, \
                              __version__), split_meta_header[0])

    def test_insert_consensus(self):
        tag = consensus_caller._DepthTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_AF:JQ_foo_DP:JQ_bar_DP:JQ_baz_DP|X:1:2:3|Y:4:5:6\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_AF:JQ_foo_DP:JQ_bar_DP:JQ_baz_DP:{0}DP_AVERAGE:{0}DP_RANGE|X:1:2:3:2.0:2.0|Y:4:5:6:5.0:2.0\n").format(consensus_caller.JQ_CONSENSUS_TAG)
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_multAlts(self):
        tag = consensus_caller._DepthTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_DP:JQ_bar_DP|0,0:1,2|0,0:3,4\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_DP:JQ_bar_DP:{0}DP_AVERAGE:{0}DP_RANGE|0,0:1,2:0.5,1.0:1.0,2.0|0,0:3,4:1.5,2.0:3.0,4.0\n").format(consensus_caller.JQ_CONSENSUS_TAG)
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_noJQ_DPTags(self):
        tag = consensus_caller._DepthTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_AF|X|Y\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_AF:{0}DP_AVERAGE:{0}DP_RANGE|X:.:.|Y:.:.\n").format(consensus_caller.JQ_CONSENSUS_TAG)
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class SomaticTagTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        split_meta_header = consensus_caller._SomaticTag().metaheader.split("\n")
        self.assertEqual('##FORMAT=<ID={0}SOM_COUNT,Number=1,Type=Integer,' \
                      'Description="Count of recognized variant callers ' \
                      'that reported confident somatic call for this '\
                      'sample-position.",Source="Jacquard",Version="{1}">'\
                      .format(consensus_caller.JQ_CONSENSUS_TAG, \
                              __version__), split_meta_header[0])

    def test_insert_consensus(self):
        tag = consensus_caller._SomaticTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP:JQ_bar_HC_SOM:JQ_baz_HC_SOM|X:2:0:1|Y:4:1:1\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP:JQ_bar_HC_SOM:JQ_baz_HC_SOM:{0}SOM_COUNT|X:2:0:1:1|Y:4:1:1:2\n").format(consensus_caller.JQ_CONSENSUS_TAG)
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_allZero(self):
        tag = consensus_caller._SomaticTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP:JQ_bar_HC_SOM:JQ_baz_HC_SOM|X:2:0:0|Y:4:0:0\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP:JQ_bar_HC_SOM:JQ_baz_HC_SOM:{0}SOM_COUNT|X:2:0:0:0|Y:4:0:0:0\n").format(consensus_caller.JQ_CONSENSUS_TAG)
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_noJQ_HC_SOMTags(self):
        tag = consensus_caller._SomaticTag()
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP|X:2|Y:4\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP:{0}SOM_COUNT|X:2:.|Y:4:.\n").format(consensus_caller.JQ_CONSENSUS_TAG)
        processedVcfRecord = VcfRecord.parse_record(line, ["SA", "SB"])
        tag.add_tag_values(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class ConsensusHelperTestCase(test_case.JacquardBaseTestCase):
    def test_add_tags(self):
        line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_VS_CALLER_REPORTED:JQ_MT_CALLER_REPORTED|0:1:1|0.2:1:1\n")
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO;JQ_CONS_SAMPLES_PASSED_COUNT=0;JQ_CONS_SAMPLES_REPORTED_COUNT=2|JQ_foo_AF:JQ_VS_CALLER_REPORTED:JQ_MT_CALLER_REPORTED:JQ_CONS_CALLERS_REPORTED_COUNT:JQ_CONS_CALLERS_PASSED_COUNT:{0}AF_AVERAGE:{0}AF_RANGE:{0}DP_AVERAGE:{0}DP_RANGE:{0}SOM_COUNT|0:1:1:2:0:0.0:.:.:.:.|0.2:1:1:2:0:0.2:.:.:.:.\n".format(consensus_caller.JQ_CONSENSUS_TAG))
        vcf_record = VcfRecord.parse_record(line, ["SA", "SB"])
        caller = consensus_caller.ConsensusCaller()
        actual = caller.add_tags(vcf_record)

        self.assertEqual(expected, actual)

    def test_get_new_metaheaders(self):
        expected = ('##FORMAT=<ID={}{},'
                    'Number=1,'
                    'Type=Integer,'
                    'Description="Count of variant callers which listed this variant in the Jacquard tagged VCF",'
                    'Source="Jacquard",'
                    'Version="{}">').format(consensus_caller.JQ_CONSENSUS_TAG,
                                            consensus_caller.JQ_REPORTED,
                                            __version__)

        caller = consensus_caller.ConsensusCaller()
        actual = caller.get_metaheaders()

        split_actual = actual[0].split("\n")
        first_meta_header = split_actual[0]

        self.assertEqual(expected, first_meta_header)
        self.assertEqual(7, len(actual))
        self.assertEqual(1, len(split_actual))

    def test_calculate_average_float(self):
        new_tags = [[0.2], [0.3], [0.5]]
        actual_tags = consensus_caller._calculate_average(new_tags)
        self.assertEquals("0.33", actual_tags)

    def test_calculate_average_int(self):
        new_tags = [[2], [3], [5]]
        actual_tags = consensus_caller._calculate_average(new_tags)
        self.assertEquals("3.33", actual_tags)

    def test_calculate_range(self):
        tag = consensus_caller._AlleleFreqTag()
        cons_freq = [[0.2], [0.4], [0.5]]
        actual_af_range = consensus_caller._calculate_range(cons_freq, tag.all_ranges)

        self.assertEquals("0.3", actual_af_range)
        self.assertEquals([0.3], tag.all_ranges)

    def test_calculate_range_oneCaller(self):
        tag = consensus_caller._AlleleFreqTag()
        cons_freq = [[0.2]]
        actual_af_range = consensus_caller._calculate_range(cons_freq, tag.all_ranges)

        self.assertEquals(".", actual_af_range)
        self.assertEquals([], tag.all_ranges)

    def test_calculate_range_multAlts(self):
        tag = consensus_caller._AlleleFreqTag()
        cons_freq = [[0.2, 0.3], [0.5, 0.7], [0.2, 0.4]]
        actual_af_range = consensus_caller._calculate_range(cons_freq, tag.all_ranges)

        self.assertEquals("0.3,0.4", actual_af_range)
        self.assertEquals([0.3, 0.4], tag.all_ranges)


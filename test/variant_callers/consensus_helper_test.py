# pylint: disable=C0103,C0301,R0904,C0111,W0212
from collections import OrderedDict
import os
import unittest

import jacquard.utils as utils
from jacquard.variant_callers import consensus_helper
from jacquard.vcf import VcfRecord

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

class MockWriter():
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
    def __init__(self, input_filepath="/foo/mockFileReader.txt", content = []):
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

class AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        split_meta_header = consensus_helper._AlleleFreqTag().metaheader.split("\n")
        self.assertEqual('##FORMAT=<ID={0}AF_AVERAGE,Number=1,Type=Float,Description="Average allele frequency across recognized variant callers that reported frequency for this position [average(JQ_*_AF)].",Source="Jacquard",Version="{1}">'.format(consensus_helper.JQ_CONSENSUS_TAG, utils.__version__), split_meta_header[0])

    def test_insert_consensus(self):
        tag = consensus_helper._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF|X:0:0.1:0.2|Y:0.2:0.3:0.4\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF:{0}AF_AVERAGE:{0}AF_RANGE|X:0:0.1:0.2:0.1:0.2|Y:0.2:0.3:0.4:0.3:0.2\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.insert_consensus(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_multAlts(self):
        tag = consensus_helper._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF|0,0:0.2,0.4|0,0:0.1,0.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF:{0}AF_AVERAGE:{0}AF_RANGE|0,0:0.2,0.4:0.1,0.2:0.2,0.4|0,0:0.1,0.3:0.05,0.15:0.1,0.3\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.insert_consensus(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_noJQ_AFTags(self):
        tag = consensus_helper._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP|X|Y\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{0}AF_AVERAGE:{0}AF_RANGE|X:.:.|Y:.:.\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.insert_consensus(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_noNullValuesInAvgCalculation(self):
        tag = consensus_helper._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF|X:0:0.1:.|Y:0.2:0.3:.\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF:{0}AF_AVERAGE:{0}AF_RANGE|X:0:0.1:.:0.05:0.1|Y:0.2:0.3:.:0.25:0.1\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.insert_consensus(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_zscore(self):
        tag = consensus_helper._AlleleFreqTag()
        content = "1\t42\t.\tA\tG\t.\tPASS\tINFO\tJQ_CONS_AF_RANGE\t0.4\t0.2"
        mock_vcf_record = MockVcfRecord(content)
        pop_mean_range = 0.5
        pop_std_range = 0.02
        tag.insert_zscore(mock_vcf_record, pop_mean_range, pop_std_range)

        self.assertTrue("JQ_CONS_AF_ZSCORE" in mock_vcf_record.format_set)
        self.assertEquals("-5.0", mock_vcf_record.sample_dict[0]["JQ_CONS_AF_ZSCORE"])
        self.assertEquals("-15.0", mock_vcf_record.sample_dict[1]["JQ_CONS_AF_ZSCORE"])

    def test_insert_zscore_nullValue(self):
        tag = consensus_helper._AlleleFreqTag()
        content = "1\t42\t.\tA\tG\t.\tPASS\tINFO\tJQ_CONS_AF_RANGE\t.\t."
        mock_vcf_record = MockVcfRecord(content)
        pop_mean_range = 0.5
        pop_std_range = 0.02
        tag.insert_zscore(mock_vcf_record, pop_mean_range, pop_std_range)

        self.assertTrue("JQ_CONS_AF_ZSCORE" in mock_vcf_record.format_set)
        self.assertEquals(".", mock_vcf_record.sample_dict[0]["JQ_CONS_AF_ZSCORE"])
        self.assertEquals(".", mock_vcf_record.sample_dict[1]["JQ_CONS_AF_ZSCORE"])

class DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        split_meta_header = consensus_helper._DepthTag().metaheader.split("\n")
        self.assertEqual('##FORMAT=<ID={0}DP_AVERAGE,Number=1,Type=Float,' \
                      'Description="Average allele frequency across ' \
                      'recognized variant callers that reported ' \
                      'frequency for this position; rounded to integer '\
                      '[round(average(JQ_*_DP))].",' \
                      'Source="Jacquard",Version="{1}">'\
                      .format(consensus_helper.JQ_CONSENSUS_TAG, \
                              utils.__version__), split_meta_header[0])

    def test_insert_consensus(self):
        tag = consensus_helper._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_AF:JQ_foo_DP:JQ_bar_DP:JQ_baz_DP|X:1:2:3|Y:4:5:6\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_AF:JQ_foo_DP:JQ_bar_DP:JQ_baz_DP:{0}DP_AVERAGE:{0}DP_RANGE|X:1:2:3:2:2|Y:4:5:6:5:2\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.insert_consensus(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_multAlts(self):
        tag = consensus_helper._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_DP:JQ_bar_DP|0,0:1,2|0,0:3,4\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_DP:JQ_bar_DP:{0}DP_AVERAGE:{0}DP_RANGE|0,0:1,2:0.5,1:1,2|0,0:3,4:1.5,2:3,4\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.insert_consensus(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_noJQ_DPTags(self):
        tag = consensus_helper._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_AF|X|Y\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_AF:{0}DP_AVERAGE:{0}DP_RANGE|X:.:.|Y:.:.\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.insert_consensus(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_zscore(self):
        tag = consensus_helper._DepthTag()
        content = "1\t42\t.\tA\tG\t.\tPASS\tINFO\tJQ_CONS_DP_RANGE\t4\t2"
        mock_vcf_record = MockVcfRecord(content)
        pop_mean_range = 5
        pop_std_range = 2
        tag.insert_zscore(mock_vcf_record, pop_mean_range, pop_std_range)

        self.assertTrue("JQ_CONS_DP_ZSCORE" in mock_vcf_record.format_set)
        self.assertEquals("-0.5", mock_vcf_record.sample_dict[0]["JQ_CONS_DP_ZSCORE"])
        self.assertEquals("-1.5", mock_vcf_record.sample_dict[1]["JQ_CONS_DP_ZSCORE"])

    def test_insert_zscore_nullValue(self):
        tag = consensus_helper._DepthTag()
        content = "1\t42\t.\tA\tG\t.\tPASS\tINFO\tJQ_CONS_DP_RANGE\t.\t."
        mock_vcf_record = MockVcfRecord(content)
        pop_mean_range = 5
        pop_std_range = 2
        tag.insert_zscore(mock_vcf_record, pop_mean_range, pop_std_range)

        self.assertTrue("JQ_CONS_DP_ZSCORE" in mock_vcf_record.format_set)
        self.assertEquals(".", mock_vcf_record.sample_dict[0]["JQ_CONS_DP_ZSCORE"])
        self.assertEquals(".", mock_vcf_record.sample_dict[1]["JQ_CONS_DP_ZSCORE"])

class SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        split_meta_header = consensus_helper._SomaticTag().metaheader.split("\n")
        self.assertEqual('##FORMAT=<ID={0}SOM_COUNT,Number=1,Type=Integer,' \
                      'Description="Count of recognized variant callers,' \
                      'which reported confident somatic call for this'\
                      'sample-position.",Source="Jacquard",Version="{1}">'\
                      .format(consensus_helper.JQ_CONSENSUS_TAG, \
                              utils.__version__), split_meta_header[0])

    def test_insert_consensus(self):
        tag = consensus_helper._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP:JQ_bar_HC_SOM:JQ_baz_HC_SOM|X:2:0:1|Y:4:1:1\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP:JQ_bar_HC_SOM:JQ_baz_HC_SOM:{0}SOM_COUNT|X:2:0:1:1|Y:4:1:1:2\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.insert_consensus(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_allZero(self):
        tag = consensus_helper._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP:JQ_bar_HC_SOM:JQ_baz_HC_SOM|X:2:0:0|Y:4:0:0\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP:JQ_bar_HC_SOM:JQ_baz_HC_SOM:{0}SOM_COUNT|X:2:0:0:0|Y:4:0:0:0\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.insert_consensus(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_insert_consensus_noJQ_HC_SOMTags(self):
        tag = consensus_helper._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP|X:2|Y:4\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_foo_DP:{0}SOM_COUNT|X:2:.|Y:4:.\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.insert_consensus(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class ConsensusHelperTestCase(unittest.TestCase):
    def test_add_tags(self):
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF:JQ_baz_AF|0:0.1:0.2|0.2:0.3:0.4\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF:JQ_baz_AF:{0}AF_AVERAGE:{0}AF_RANGE:{0}DP_AVERAGE:{0}DP_RANGE|0:0.1:0.2:0.1:0.2:.:.|0.2:0.3:0.4:0.3:0.2:.:.\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        vcf_record = VcfRecord(line)
        cons_help = consensus_helper.ConsensusHelper()
        actual = cons_help.add_tags(vcf_record)

        self.assertEqual(expected, actual)

    def test_get_new_metaheaders(self):
        expected = ('##FORMAT=<ID={0}AF_AVERAGE,Number=1,Type=Float,'+
        'Description="Average allele frequency across recognized variant '+
        'callers that reported frequency for this position '+
        '[average(JQ_*_AF)].",Source="Jacquard",'+
        'Version="{1}">').format(consensus_helper.JQ_CONSENSUS_TAG,
                                 utils.__version__)

        cons_help = consensus_helper.ConsensusHelper()
        actual = cons_help.get_new_metaheaders()

        split_actual = actual[0].split("\n")
        first_meta_header = split_actual[0]

        self.assertEqual(expected, first_meta_header)
        self.assertEqual(2, len(actual))
        self.assertEqual(3, len(split_actual))

    def test_calculate_average_float(self):
        desired_tags = [[0.2],[0.3],[0.5]]
        actual_tags = consensus_helper._calculate_average(desired_tags)
        self.assertEquals("0.33", actual_tags)

    def test_calculate_average_int(self):
        desired_tags = [[2],[3],[5]]
        actual_tags = consensus_helper._calculate_average(desired_tags)
        self.assertEquals("3.33", actual_tags)

    def test_calculate_range(self):
        tag = consensus_helper._AlleleFreqTag()
        cons_freq = [[0.2], [0.4], [0.5]]
        actual_af_range = consensus_helper._calculate_range(cons_freq, tag.all_ranges)

        self.assertEquals("0.3", actual_af_range)
        self.assertEquals([0.3], tag.all_ranges)

    def test_calculate_range_oneCaller(self):
        tag = consensus_helper._AlleleFreqTag()
        cons_freq = [[0.2]]
        actual_af_range = consensus_helper._calculate_range(cons_freq, tag.all_ranges)

        self.assertEquals(".", actual_af_range)
        self.assertEquals([], tag.all_ranges)

    def test_calculate_range_multAlts(self):
        tag = consensus_helper._AlleleFreqTag()
        cons_freq = [[0.2, 0.3], [0.5, 0.7], [0.2, 0.4]]
        actual_af_range = consensus_helper._calculate_range(cons_freq, tag.all_ranges)

        self.assertEquals("0.3,0.4", actual_af_range)
        self.assertEquals([0.3, 0.4], tag.all_ranges)

    def test_get_pop_values(self):
        all_ranges = {0: [0.2, 0.6, 0.8, 0.1, 0.0, 0.3]}
        pop_mean_range, pop_std_range = consensus_helper._calculate_population_values(all_ranges)

        self.assertEquals(0.33, pop_mean_range)
        self.assertEquals(0.28, pop_std_range)


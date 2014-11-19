# pylint: disable=C0103,C0301,R0904,C0111
import os
import unittest

import jacquard.utils as utils
from jacquard.variant_callers import consensus_helper
from jacquard.vcf import VcfRecord

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
        self.assertEqual('##FORMAT=<ID={0}AF_AVERAGE,Number=1,Type=Float,Description="Average allele frequency across recognized variant callers that reported frequency for this position [average(JQ_*_AF)].", Source="Jacquard", Version="{1}">'.format(consensus_helper.JQ_CONSENSUS_TAG, utils.__version__), consensus_helper._AlleleFreqTag().metaheader)

    def test_format(self):
        tag = consensus_helper._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF|X:0:0.1:0.2|Y:0.2:0.3:0.4\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF:{0}AF_AVERAGE|X:0:0.1:0.2:0.1|Y:0.2:0.3:0.4:0.3\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_multAlts(self):
        tag = consensus_helper._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF|0,0:0.2,0.4|0,0:0.1,0.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF:{0}AF_AVERAGE|0,0:0.2,0.4:0.1,0.2|0,0:0.1,0.3:0.05,0.15\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_unequalMultAlts(self):
        tag = consensus_helper._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF|0,0:0.4|0,0:0.1\n".replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        self.assertRaisesRegexp(utils.JQException, "Inconsistent number of mult-alts found in VCF file", tag.format, processedVcfRecord)

    def test_format_noJQ_AFTags(self):
        tag = consensus_helper._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP|X|Y\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:{0}AF_AVERAGE|X:.|Y:.\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)

        self.assertEquals(expected, processedVcfRecord.asText())

class ConsensusHelperTestCase(unittest.TestCase):
    def test_add_tags(self):
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF|X:0:0.1:0.2|Y:0.2:0.3:0.4\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_DP:JQ_foo_AF:JQ_bar_AF:JQ_baz_AF:{0}AF_AVERAGE|X:0:0.1:0.2:0.1|Y:0.2:0.3:0.4:0.3\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        vcf_record = VcfRecord(line)
        cons_help = consensus_helper.ConsensusHelper()
        actual = cons_help.add_tags(vcf_record)
        self.assertEqual(expected, actual)
    
    def test_get_new_metaheaders(self):
        expected = [('##FORMAT=<ID={0}AF_AVERAGE,Number=1,Type=Float,'+
        'Description="Average allele frequency across recognized variant '+
        'callers that reported frequency for this position '+
        '[average(JQ_*_AF)].", Source="Jacquard", '+
        'Version="{1}">').format(consensus_helper.JQ_CONSENSUS_TAG, 
                                 utils.__version__)]

        cons_help = consensus_helper.ConsensusHelper()
        actual = cons_help.get_new_metaheaders()
        self.assertEqual(expected, actual)

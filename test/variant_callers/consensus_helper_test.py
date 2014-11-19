'''
Created on Nov 18, 2014

@author: kmeng
'''
# pylint: disable=C0103,C0301,R0904
from collections import OrderedDict,defaultdict
import os
import unittest

from jacquard.utils import __version__,jq_af_tag,jq_dp_tag,jq_somatic_tag,\
    JQException
import jacquard.utils as utils
from jacquard.variant_callers import consensus_helper

from jacquard.vcf import VcfRecord, FileReader
from testfixtures import TempDirectory

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


class ConsensusHelperTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}AF_AVERAGE,Number=1,Type=Float,Description="Average allele frequency across recognized variant callers that reported frequency for this position [average(JQ_*_AF)].", Source="Jacquard", Version="{1}">'.format(consensus_helper.JQ_CONSENSUS_TAG, utils.__version__), consensus_helper._AlleleFreqTag().metaheader)

    def test_format(self):
        tag = consensus_helper._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF:JQ_baz_AF|0:1:2|2:3:4\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF:JQ_baz_AF:{0}AF_AVERAGE|0:1:2:1|2:3:4:3\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        
    def test_format_multAlts(self):
        tag = consensus_helper._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF|0,0:2,4|0,0:1,3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|JQ_foo_AF:JQ_bar_AF:{0}AF_AVERAGE|0,0:2,4:1,2|0,0:1,3:0.5,1.5\n".format(consensus_helper.JQ_CONSENSUS_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())
        

        
        
#!/usr/bin/python2.7
from collections import OrderedDict
import unittest

from bin.consensus import add_consensus, process_line, calculate_consensus, create_consensus_dict, get_consensus_som, get_consensus

class ConsensusTestCase(unittest.TestCase):
    def test_addConsensusSomatic(self):
        reader = MockReader("1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_AF_VS\t1:0:0.0\t0:1:1.2\t1:3:1.0\n")
        writer = MockWriter()
        output_file = "output"
        add_consensus(reader, writer, output_file)
        self.assertEquals(['##FORMAT=<ID=JQ_SOM_SUM,Number=1,Type=Integer,Description="Jacquard consensus somatic call = sum(JQ_SOM_MT, JQ_SOM_SK, JQ_SOM_VS)">',
                            '##FORMAT=<ID=JQ_AF_AVERAGE,Number=A,Type=Integer,Description="Jacquard consensus somatic call = average(JQ_AF_MT, JQ_AF_SK, JQ_AF_VS)">',
                            '##FORMAT=<ID=JQ_DP_AVERAGE,Number=1,Type=Integer,Description="Jacquard consensus depth = average(JQ_DP_MT, JQ_DP_SK, JQ_DP_VS)">',
                            '1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_AF_VS:JQ_AF_AVERAGE:JQ_SOM_SUM:JQ_DP_AVERAGE\t1:0:0.0:0.0:1:0\t0:1:1.2:1.2:0:0\t1:3:1.0:1.0:1:0'], writer.lines())
        
    def test_addConsensusSomatic_multipleCallers(self):
        reader = MockReader("1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_SOM_MT\t1:0:1\t0:1:0\t1:3:0\n")
        writer = MockWriter()
        output_file = "output"
        add_consensus(reader, writer, output_file)
        self.assertEquals(['##FORMAT=<ID=JQ_SOM_SUM,Number=1,Type=Integer,Description="Jacquard consensus somatic call = sum(JQ_SOM_MT, JQ_SOM_SK, JQ_SOM_VS)">',
                           '##FORMAT=<ID=JQ_AF_AVERAGE,Number=A,Type=Integer,Description="Jacquard consensus somatic call = average(JQ_AF_MT, JQ_AF_SK, JQ_AF_VS)">',
                           '##FORMAT=<ID=JQ_DP_AVERAGE,Number=1,Type=Integer,Description="Jacquard consensus depth = average(JQ_DP_MT, JQ_DP_SK, JQ_DP_VS)">',
                           '1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_SOM_MT:JQ_AF_AVERAGE:JQ_SOM_SUM:JQ_DP_AVERAGE\t1:0:1:0:2:0\t0:1:0:0:0:0\t1:3:0:0:1:0'], writer.lines())
        
    def test_calculateConsensus(self):
        combined_dict = OrderedDict([("DP", "34"), ("JQ_SOM_VS", "1"), ("JQ_SOM_MT", "0"), ("JQ_AF_VS", "0.5"), ("JQ_DP_VS", "53")])
        actual_dict = calculate_consensus(combined_dict)
        expected_dict = OrderedDict([("DP", "34"), ("JQ_SOM_VS", "1"), ("JQ_SOM_MT", "0"), ("JQ_AF_VS", "0.5"), ("JQ_DP_VS", "53"), ("JQ_AF_AVERAGE", "0.5"), ("JQ_SOM_SUM", "1"), ("JQ_DP_AVERAGE", "53.0")])
        self.assertEquals(expected_dict, actual_dict)
        
    def test_calculateConsensus_multipleCallers(self):
        combined_dict = OrderedDict([("DP", "34"), ("JQ_SOM_VS", "1"), ("JQ_SOM_MT", "1"), ("JQ_SOM_FOO", "1"), ("JQ_AF_VS", "0.5"), ("JQ_AF_MT", "0.43"), ("JQ_DP_MT", "32"), ("JQ_DP_VS", "23")])
        actual_dict = calculate_consensus(combined_dict)
        expected_dict = OrderedDict([("DP", "34"), ("JQ_SOM_VS", "1"), ("JQ_SOM_MT", "1"), ("JQ_SOM_FOO", "1"), ("JQ_AF_VS", "0.5"), ("JQ_AF_MT", "0.43"), ("JQ_DP_MT", "32"), ("JQ_DP_VS", "23"), ("JQ_AF_AVERAGE", "0"), ("JQ_AF_AVERAGE", "0.465"), ("JQ_SOM_SUM", "3"), ("JQ_DP_AVERAGE", "27.5")])
        self.assertEquals(expected_dict, actual_dict)
        
    def test_calculateConsensus_multipleAlts(self):
        combined_dict = OrderedDict([("DP", "34"), ("JQ_SOM_VS", "1,0"), ("JQ_SOM_MT", "1,1"), ("JQ_SOM_FOO", "1,0"), ("JQ_AF_VS", "0.5,0.0"), ("JQ_AF_MT", "0.43,0.3"), ("JQ_DP_MT", "32,45"), ("JQ_DP_VS", "23,10")])
        actual_dict = calculate_consensus(combined_dict)
        expected_dict = OrderedDict([("DP", "34"), ("JQ_SOM_VS", "1,0"), ("JQ_SOM_MT", "1,1"), ("JQ_SOM_FOO", "1,0"), ("JQ_AF_VS", "0.5,0.0"), ("JQ_AF_MT", "0.43,0.3"), ("JQ_DP_MT", "32,45"), ("JQ_DP_VS", "23,10"), ("JQ_AF_AVERAGE", "0.465,0.15"), ("JQ_SOM_SUM", "3,1"), ("JQ_DP_AVERAGE", "27.5,27.5")])
        self.assertEquals(expected_dict, actual_dict)
        
    def test_processLine(self):
        line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_SOM_MT:JQ_AF_VS:JQ_DP_VS\t1:0:1:0.2:43\t0:1:0:0.42:0\t1:3:0:0.0:12\n"
        new_line = process_line(line)
        expected_line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_SOM_MT:JQ_AF_VS:JQ_DP_VS:JQ_AF_AVERAGE:JQ_SOM_SUM:JQ_DP_AVERAGE\t1:0:1:0.2:43:0.2:2:43.0\t0:1:0:0.42:0:0.42:0:0.0\t1:3:0:0.0:12:0.0:1:12.0\n"
        self.assertEquals(expected_line, new_line)
    
    def test_processLine_noJQTags(self):
        line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tFOO:DP\t1:0\t0:2\t1:3\n"
        new_line = process_line(line)
        expected_line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tFOO:DP:JQ_AF_AVERAGE:JQ_SOM_SUM:JQ_DP_AVERAGE\t1:0:0:0:0\t0:2:0:0:0\t1:3:0:0:0\n"
        self.assertEquals(expected_line, new_line)
    
    def test_createConsensusDict_som(self):
        input_dict = OrderedDict([("JQ_SOM_VS", "0,1"), ("DP", "23")])
        consensus_dict = create_consensus_dict("JQ_SOM_VS", "0,1", input_dict, {}, "int")
        expected_dict = {0: 0, 1: 1}
        self.assertEquals(expected_dict, consensus_dict)
        
    def test_createConsensusDict_af(self):
        input_dict = OrderedDict([("JQ_AF_VS", "0.02,0.23"), ("DP", "23")])
        consensus_dict = create_consensus_dict("JQ_AF_VS", "0.02,0.23", input_dict, {}, "float")
        expected_dict = {0: 0.02, 1: 0.23}
        self.assertEquals(expected_dict, consensus_dict)
        
    def test_createConsensusDict_dp(self):
        input_dict = OrderedDict([("JQ_DP_VS", "2,23"), ("DP", "23")])
        consensus_dict = create_consensus_dict("JQ_DP_VS", "2,23", input_dict, {}, "float")
        expected_dict = {0: 2.0, 1: 23.0}
        self.assertEquals(expected_dict, consensus_dict)
        
    def test_getConsensusSom(self):
        field_dict = {0:1, 1:4, 2:0}
        consensus = get_consensus_som(field_dict)
        self.assertEquals("1,4,0", consensus)
        
        field_dict = {0:1}
        consensus = get_consensus_som(field_dict)
        self.assertEquals("1", consensus)
    
    def test_getConsensus(self):
        tags = ["1", "2"]
        field_dict = {0:1.0, 1:4.4, 2:0.1}
        consensus = get_consensus(tags, field_dict)
        self.assertEquals("0.5,2.2,0.05", consensus)
        
        tags = ["1"]
        field_dict = {0:1.0}
        consensus = get_consensus(tags, field_dict)
        self.assertEquals("1.0", consensus)
        
class MockWriter():
    def __init__(self):
        self._content = []
        self.wasClosed = False

    def write(self, content):
        self._content.extend(content.splitlines())
        
    def lines(self):
        return self._content

    def close(self):
        self.wasClosed = True
        
class MockReader():
    def __init__(self, content):
        lines = [line + "\n" for line in content.split("\n") if line != ""]
        self._iter = lines.__iter__()
        self.wasClosed = False

    def __iter__(self):
        return self._iter

    def close(self):
        self.wasClosed=True
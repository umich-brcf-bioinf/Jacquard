from collections import OrderedDict
import numpy
import unittest

from jacquard.consensus import iterate_file, add_consensus, process_line, calculate_consensus, create_consensus_dict, get_consensus_som, get_consensus, add_zscore, calculate_zscore

class ConsensusTestCase(unittest.TestCase):
    def test_iterateFile_consensus(self):
        reader = MockReader("##FOO\n#CHROM\n1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_HC_SOM_VS:DP:JQ_AF_VS:JQ_AF_AVERAGE:JQ_DP_AVERAGE\t1:0:0.0:0.0:0\t0:1:1.2:1.2:0\t1:3:1.0:1.0:0\n")
        writer = MockWriter()
        output_file = "output"
        af_range = []
        dp_range = []
        type = "consensus"
        meta_headers, header, lines = iterate_file(reader, writer, output_file, af_range, dp_range, type)
        
        self.assertEquals(["##FOO\n"], meta_headers)
        self.assertEquals("#CHROM\n", header)
        self.assertEquals(["1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_HC_SOM_VS:DP:JQ_AF_VS:JQ_AF_AVERAGE:JQ_DP_AVERAGE:JQ_SOM_SUM:JQ_AF_RANGE:JQ_DP_RANGE\t1:0:0.0:0.0:0.0:1:0.0:0\t0:1:1.2:1.2:0.0:0:0.0:0\t1:3:1.0:1.0:0.0:1:0.0:0\n"], lines)
        
    def test_iterateFile_zscore(self):
        reader = MockReader("##FOO\n#CHROM\n1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_HC_SOM_VS:DP:JQ_AF_VS:JQ_AF_AVERAGE:JQ_DP_AVERAGE:JQ_SOM_SUM:JQ_AF_RANGE:JQ_DP_RANGE\t1:0:0.0:0.0:0.0:1:0.0:0\t0:1:1.2:1.2:0.0:0:0.0:0\t1:3:1.0:1.0:0.0:1:0.0:0\n")
        writer = MockWriter()
        output_file = "output"
        af_range = [0.2, 0.45]
        dp_range = [45.1, 86.0]
        type = "zscore"
        meta_headers, header, lines = iterate_file(reader, writer, output_file, af_range, dp_range, type)
        
        self.assertEquals(["##FOO\n"], meta_headers)
        self.assertEquals("#CHROM\n", header)
        self.assertEquals(["1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_HC_SOM_VS:DP:JQ_AF_VS:JQ_AF_AVERAGE:JQ_DP_AVERAGE:JQ_SOM_SUM:JQ_AF_RANGE:JQ_DP_RANGE:JQ_AF_RANGE_ZSCORE:JQ_DP_RANGE_ZSCORE\t1:0:0.0:0.0:0.0:1:0.0:0:-2.6:-3.21\t0:1:1.2:1.2:0.0:0:0.0:0:-2.6:-3.21\t1:3:1.0:1.0:0.0:1:0.0:0:-2.6:-3.21\n"], lines)
        
    def test_addConsensusSomatic(self):
        meta_headers = []
        header = "#CHROM"
        lines = ["1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_AF_VS:JQ_AF_AVERAGE:JQ_SOM_SUM:JQ_DP_AVERAGE:JQ_AF_RANGE:JQ_DP_RANGE\t1:0:0.0:0.0:1:0:0:0\t0:1:1.2:1.2:0:0:0:0\t1:3:1.0:1.0:1:0:0:0"]
        writer = MockWriter()
        output_file = "output"
        add_consensus(meta_headers, header, lines, writer, output_file)
        self.assertEquals(['##FORMAT=<ID=JQ_SOM_SUM,Number=1,Type=Integer,Description="Jacquard consensus somatic call = sum(JQ_HC_SOM_*)">',
                            '##FORMAT=<ID=JQ_AF_AVERAGE,Number=A,Type=Integer,Description="Jacquard consensus somatic call = average(JQ_AF_*)">',
                            '##FORMAT=<ID=JQ_DP_AVERAGE,Number=1,Type=Integer,Description="Jacquard consensus depth = average(JQ_DP_*)">',
                            '#CHROM',
                            '1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_AF_VS:JQ_AF_AVERAGE:JQ_SOM_SUM:JQ_DP_AVERAGE:JQ_AF_RANGE:JQ_DP_RANGE\t1:0:0.0:0.0:1:0:0:0\t0:1:1.2:1.2:0:0:0:0\t1:3:1.0:1.0:1:0:0:0'], writer.lines())
        
    def test_addConsensusSomatic_multipleCallers(self):
        meta_headers = []
        header = "#CHROM"
        lines = ["1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_SOM_MT:JQ_AF_AVERAGE:JQ_SOM_SUM:JQ_DP_AVERAGE:JQ_AF_RANGE:JQ_DP_RANGE\t1:0:1:0:2:0:0:0\t0:1:0:0:0:0:0:0\t1:3:0:0:1:0:0:0\n"]
        writer = MockWriter()
        output_file = "output"
        add_consensus(meta_headers, header, lines, writer, output_file)
        self.assertEquals(['##FORMAT=<ID=JQ_SOM_SUM,Number=1,Type=Integer,Description="Jacquard consensus somatic call = sum(JQ_HC_SOM_*)">',
                           '##FORMAT=<ID=JQ_AF_AVERAGE,Number=A,Type=Integer,Description="Jacquard consensus somatic call = average(JQ_AF_*)">',
                           '##FORMAT=<ID=JQ_DP_AVERAGE,Number=1,Type=Integer,Description="Jacquard consensus depth = average(JQ_DP_*)">',
                           '#CHROM', 
                           '1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:DP:JQ_SOM_MT:JQ_AF_AVERAGE:JQ_SOM_SUM:JQ_DP_AVERAGE:JQ_AF_RANGE:JQ_DP_RANGE\t1:0:1:0:2:0:0:0\t0:1:0:0:0:0:0:0\t1:3:0:0:1:0:0:0'], writer.lines())
        
    def test_calculateConsensus(self):
        combined_dict = OrderedDict([("DP", "34"), ("JQ_HC_SOM_VS", "1"), ("JQ_HC_SOM_MT", "0"), ("JQ_AF_VS", "0.5"), ("JQ_DP_VS", "53")])
        actual_dict, af, dp = calculate_consensus(combined_dict, [], [])
        expected_dict = OrderedDict([("DP", "34"), ("JQ_HC_SOM_VS", "1"), ("JQ_HC_SOM_MT", "0"), ("JQ_AF_VS", "0.5"), ("JQ_DP_VS", "53"), ("JQ_AF_AVERAGE", "0.5"), ("JQ_SOM_SUM", "1"), ("JQ_DP_AVERAGE", "53.0"), ('JQ_AF_RANGE', '0'), ('JQ_DP_RANGE', '0')])
        self.assertEquals(expected_dict, actual_dict)
        self.assertEquals([0], af)
        self.assertEquals([0], dp)
        
    def test_calculateConsensus_multipleCallers(self):
        combined_dict = OrderedDict([("DP", "34"), ("JQ_HC_SOM_VS", "1"), ("JQ_HC_SOM_MT", "1"), ("JQ_HC_SOM_FOO", "1"), ("JQ_AF_VS", "0.5"), ("JQ_AF_MT", "0.43"), ("JQ_DP_MT", "32"), ("JQ_DP_VS", "23")])
        actual_dict, af, dp = calculate_consensus(combined_dict, [], [])
        expected_dict = OrderedDict([("DP", "34"), ("JQ_HC_SOM_VS", "1"), ("JQ_HC_SOM_MT", "1"), ("JQ_HC_SOM_FOO", "1"), ("JQ_AF_VS", "0.5"), ("JQ_AF_MT", "0.43"), ("JQ_DP_MT", "32"), ("JQ_DP_VS", "23"), ("JQ_AF_AVERAGE", "0"), ("JQ_AF_AVERAGE", "0.47"), ("JQ_SOM_SUM", "3"), ("JQ_DP_AVERAGE", "27.5"), ('JQ_AF_RANGE', '0.07'), ('JQ_DP_RANGE', '9.0')])
        self.assertEquals(expected_dict, actual_dict)
        self.assertEquals([0.07], af)
        self.assertEquals([9.0], dp)
        
    def test_processLine_consensus(self):
        line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_HC_SOM_VS:DP:JQ_HC_SOM_MT:JQ_AF_VS:JQ_DP_VS\t1:0:1:0.2:43\t0:1:0:0.42:0\t1:3:0:0.0:12\n"
        new_line = process_line(line, [], "", "", [], "", "", "consensus")
        expected_line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_HC_SOM_VS:DP:JQ_HC_SOM_MT:JQ_AF_VS:JQ_DP_VS:JQ_AF_AVERAGE:JQ_SOM_SUM:JQ_DP_AVERAGE:JQ_AF_RANGE:JQ_DP_RANGE\t1:0:1:0.2:43:0.2:2:43.0:0:0\t0:1:0:0.42:0:0.42:0:0.0:0:0\t1:3:0:0.0:12:0.0:1:12.0:0:0\n"
        self.assertEquals(expected_line, new_line)
    
    def test_processLine_consensusNoJQTags(self):
        line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tFOO:DP\t1:0\t0:2\t1:3\n"
        new_line = process_line(line, [], "", "", [], "", "", "consensus")
        expected_line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tFOO:DP:JQ_AF_AVERAGE:JQ_SOM_SUM:JQ_DP_AVERAGE:JQ_AF_RANGE:JQ_DP_RANGE\t1:0:0:0:0:0:0\t0:2:0:0:0:0:0\t1:3:0:0:0:0:0\n"
        self.assertEquals(expected_line, new_line)
    
    def test_processLine_zscore(self):
        af_range = [0.1, 1.2, 1.0]
        af_mean = sum(af_range)/len(af_range) if af_range != [] else ""
        af_std = numpy.std(af_range) if af_range != [] else ""
        
        dp_range = [10.0, 1.0, 3.0]
        dp_mean = sum(dp_range)/len(dp_range) if dp_range != [] else ""
        dp_std = numpy.std(dp_range) if dp_range != [] else ""
        
        line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:JQ_DP_VS:JQ_AF_VS:JQ_AF_RANGE:JQ_DP_RANGE\t1:10.0:0.1:0.1:10.0\t0:1.0:1.2:1.2:1.0\t1:3.0:1.0:1.0:3.0\n"
        actual_line = process_line(line, af_range, af_mean, af_std, dp_range, dp_mean, dp_std, "zscore")
        
        expected_line = "1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:JQ_DP_VS:JQ_AF_VS:JQ_AF_RANGE:JQ_DP_RANGE:JQ_AF_RANGE_ZSCORE:JQ_DP_RANGE_ZSCORE\t1:10.0:0.1:0.1:10.0:-1.39:1.38\t0:1.0:1.2:1.2:1.0:0.91:-0.95\t1:3.0:1.0:1.0:3.0:0.49:-0.43\n"
        self.assertEquals(expected_line, actual_line)
    
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
        
    def test_addZscore(self):
        meta_headers = []
        header = "#CHROM"
        lines = ["1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:JQ_DP_VS:JQ_AF_VS:JQ_AF_RANGE:JQ_DP_RANGE:JQ_AF_RANGE_ZSCORE:JQ_DP_RANGE_ZSCORE\t1:10.0:0.1:0.1:10.0:-1.39:1.38\t0:1.0:1.2:1.2:1.0:0.91:-0.95\t1:3.0:1.0:1.0:3.0:0.49:-0.43\n"]
        writer = MockWriter()
        output_file = "output"
        af_range = [0.1, 1.2, 1.0]
        dp_range = [10.0, 1.0, 3.0]
        add_zscore(meta_headers, header, lines, writer, output_file, af_range, dp_range)
        self.assertEquals(['##FORMAT=<ID=JQ_AF_RANGE_ZSCORE,Number=A,Type=Integer,Description="Jacquard measure of consistency of allele frequencies among callers = (sample AF range - population mean AF range)/standard dev(population AF range)">',
                           '##jacquard.consensus.JQ_AF_RANGE_ZSCORE.mean_AF_range=0.77',
                            '##jacquard.consensus.JQ_AF_RANGE_ZSCORE.standard_deviation=0.48', 
                            '##FORMAT=<ID=JQ_DP_RANGE_ZSCORE,Number=A,Type=Integer,Description="Jacquard measure of consistency of depth among callers = (sample DP range - population mean DP range)/standard dev(population DP range)">', 
                            '##jacquard.consensus.JQ_DP_RANGE_ZSCORE.mean_DP_range=4.67', 
                            '##jacquard.consensus.JQ_DP_RANGE_ZSCORE.standard deviation_DP_range=3.86', 
                            '#CHROM', 
                            '1\t1344\t.\tA\tT\t.\t.\tfoo\tJQ_SOM_VS:JQ_DP_VS:JQ_AF_VS:JQ_AF_RANGE:JQ_DP_RANGE:JQ_AF_RANGE_ZSCORE:JQ_DP_RANGE_ZSCORE\t1:10.0:0.1:0.1:10.0:-1.39:1.38\t0:1.0:1.2:1.2:1.0:0.91:-0.95\t1:3.0:1.0:1.0:3.0:0.49:-0.43'], writer.lines())
        
    def test_calculateZscore(self):
        af_range = [0.1, 1.2, 1.0]
        af_mean = sum(af_range)/len(af_range) if af_range != [] else ""
        af_std = numpy.std(af_range) if af_range != [] else ""
        
        dp_range = [10.0, 1.0, 3.0]
        dp_mean = sum(dp_range)/len(dp_range) if dp_range != [] else ""
        dp_std = numpy.std(dp_range) if dp_range != [] else ""
        
        combined_dict = {"JQ_AF_RANGE" : "0.1", "JQ_DP_RANGE" : "10.0"}
        combined_dict = calculate_zscore(af_mean, af_std, dp_mean, dp_std, combined_dict)
        
        self.assertEquals({'JQ_DP_RANGE_ZSCORE': '1.38', 'JQ_AF_RANGE_ZSCORE': '-1.39', 'JQ_AF_RANGE': '0.1', 'JQ_DP_RANGE': '10.0'}, combined_dict)
        
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
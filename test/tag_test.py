#!/usr/bin/python2.7
from collections import OrderedDict
import glob
import os
from os import listdir
from StringIO import StringIO
import sys
from stat import *
import testfixtures
from testfixtures import TempDirectory
import unittest
from bin.tag import Varscan, Mutect, Unknown, Varscan_AlleleFreqTag, Varscan_DepthTag, Varscan_SomaticTag, Mutect_AlleleFreqTag, Mutect_DepthTag, Mutect_SomaticTag, Strelka_AlleleFreqTag, Strelka_DepthTag, Strelka_SomaticTag, LineProcessor, FileProcessor, tag_files, determine_file_types, print_file_types
from bin.jacquard_utils import __version__

class Varscan_AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_VS,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ),Source="Jacquard",Version={0}>\n'.format(__version__), Varscan_AlleleFreqTag().metaheader)
                
    def test_format_missingAFTag(self):
        tag = Varscan_AlleleFreqTag()
        info_string = ""
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("alt", "filter", info_string, format_dict, 0))
                
    def test_format_rounds(self):
        tag = Varscan_AlleleFreqTag()
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20%'), ('JQ_AF_VS', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.0%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.0%'), ('JQ_AF_VS', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.4%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.4%'), ('JQ_AF_VS', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.5%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.5%'), ('JQ_AF_VS', '0.21')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.6%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.6%'), ('JQ_AF_VS', '0.21')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:100%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '100%'), ('JQ_AF_VS', '1.0')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:100.0%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '100.0%'), ('JQ_AF_VS', '1.0')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.4%,38.075%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.4%,38.075%'), ('JQ_AF_VS', '0.2,0.38')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FREQ".split(":"), "1:20.4%,38.075%,27.843%".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FREQ', '20.4%,38.075%,27.843%'), ('JQ_AF_VS', '0.2,0.38,0.28')]), tag.format("alt", "filter", "", format_dict, 0))
        
class Varscan_DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_DP_VS,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP),Source="Jacquard",Version={0}>\n'.format(__version__), Varscan_DepthTag().metaheader)
                
    def test_format_missingDPTag(self):
        tag = Varscan_DepthTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("alt", "filter", "", format_dict, 0))
                
    def test_format(self):
        tag = Varscan_DepthTag()
        format_dict = OrderedDict(zip("A:DP".split(":"), "1:42".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('DP', '42'), ('JQ_DP_VS', '42')]), tag.format("alt", "filter", "", format_dict, 0))

class Varscan_SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_SOM_VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1= somatic (based on SOMATIC info tag and if sample is TUMOR),Source="Jacquard",Version={0}>\n'.format(__version__), Varscan_SomaticTag().metaheader)
                
    def test_format_missingSSInfoTag(self):
        tag = Varscan_SomaticTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(["A", "B", "JQ_HC_SOM_VS"], tag.format("alt", "filter", "INFO", format_dict, 0).keys())
        self.assertEqual(["1", "2", "0"], tag.format("alt", "filter", "INFO", format_dict, 0).values())
                
    def test_format(self):
        tag = Varscan_SomaticTag()
        format_dict = OrderedDict([("A", "1")])
        self.assertEqual(OrderedDict([("A", "1"), ("JQ_HC_SOM_VS", "0")]), tag.format("alt", "filter", "INFO;SS=2", format_dict, 0))
        self.assertEqual(OrderedDict([("A", "1"), ("JQ_HC_SOM_VS", "0")]), tag.format("alt", "filter", "INFO;SS=2", format_dict, 1))
        self.assertEqual(OrderedDict([("A", "1"), ("JQ_HC_SOM_VS", "1")]), tag.format("alt", "filter", "INFO;SS=2;JQ_HC_VS", format_dict, 1))
        
        format_dict = OrderedDict([("A", "1")])
        self.assertEqual(OrderedDict([("A", "1",), ("JQ_HC_SOM_VS", "0")]), tag.format("alt", "filter", "INFO", format_dict, 0))
        self.assertEqual(OrderedDict([("A", "1"), ("JQ_HC_SOM_VS", "0")]), tag.format("alt", "filter", "INFO", format_dict, 1))

class Mutect_AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_MT,Number=A,Type=Float,Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA),Source="Jacquard",Version={0}>\n'.format(__version__), Mutect_AlleleFreqTag().metaheader)
                
    def test_format_missingAFTag(self):
        tag = Mutect_AlleleFreqTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("alt", "filter", "", format_dict, 0))
                
    def test_format_rounds(self):
        tag = Mutect_AlleleFreqTag()
        format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.2".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.2'), ('JQ_AF_MT', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.20".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.20'), ('JQ_AF_MT', '0.20')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.204".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.204'), ('JQ_AF_MT', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.205".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.205'), ('JQ_AF_MT', '0.21')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.206".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.206'), ('JQ_AF_MT', '0.21')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FA".split(":"), "1:1.0".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FA', '1.0'), ('JQ_AF_MT', '1.0')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:FA".split(":"), "1:1.00".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FA', '1.00'), ('JQ_AF_MT', '1.00')]), tag.format("alt", "filter", "", format_dict, 0))

        format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.204,0.3807".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.204,0.3807'), ('JQ_AF_MT', '0.2,0.38')]), tag.format("alt", "filter", "", format_dict, 0))
 
        format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.204,0.3807,0.2784".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.204,0.3807,0.2784'), ('JQ_AF_MT', '0.2,0.38,0.28')]), tag.format("alt", "filter", "", format_dict, 0))

class Mutect_DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_DP_MT,Number=1,Type=Float,Description="Jacquard depth for MuTect (based on DP),Source="Jacquard",Version={0}>\n'.format(__version__), Mutect_DepthTag().metaheader)
                
    def test_format_missingDPTag(self):
        tag = Mutect_DepthTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("alt", "filter", "", format_dict, 0))
                
    def test_format(self):
        tag = Mutect_DepthTag()
        format_dict = OrderedDict(zip("A:DP".split(":"), "1:42".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('DP', '42'), ('JQ_DP_MT', '42')]), tag.format("alt", "filter", "", format_dict, 0))

class Mutect_SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_SOM_MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag),Source="Jacquard",Version={0}>\n'.format(__version__), Mutect_SomaticTag().metaheader)
                
    def test_format_missingSSTag(self):
        tag = Mutect_SomaticTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2'), ("JQ_HC_SOM_MT", "0")]), tag.format("alt", "filter", "", format_dict, 0))
                
    def test_format(self):
        tag = Mutect_SomaticTag()
        format_dict = OrderedDict(zip( "A:SS".split(":"), "1:2".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('SS', '2'), ('JQ_HC_SOM_MT', '1')]), tag.format("alt", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip( "A:SS".split(":"), "1:1".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('SS', '1'), ('JQ_HC_SOM_MT', '0')]), tag.format("alt", "filter", "", format_dict, 0))

class Strelka_AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_MT,Number=A,Type=Float,Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA),Source="Jacquard",Version={0}>\n'.format(__version__), Mutect_AlleleFreqTag().metaheader)
                
    def test_format_missingAFTag(self):
        tag = Strelka_AlleleFreqTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("C", "filter", "", format_dict, 0))
                
    def test_format_rounds(self):
        tag = Strelka_AlleleFreqTag()
        format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:0.2,0.5:0.2,0.5:0.2,0.5:0.2,0.5:".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('AU', '0.2,0.5'), ('CU', '0.2,0.5'), ('GU', '0.2,0.5'), ('TU', '0.2,0.5'), ('JQ_AF_SK', '0.25')]), tag.format("C", "filter", "", format_dict, 0))
#         
        format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:0.2,0.6:0.2,0.6:0.2,0.5:0.2,0.5:".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('AU', '0.2,0.6'), ('CU', '0.2,0.6'), ('GU', '0.2,0.5'), ('TU', '0.2,0.5'), ('JQ_AF_SK', '0.27')]), tag.format("C", "filter", "", format_dict, 0))
#         
        format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:0.2,0.6:0.2,0.8:0.2,0.5:0.2,0.5:".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('AU', '0.2,0.6'), ('CU', '0.2,0.8'), ('GU', '0.2,0.5'), ('TU', '0.2,0.5'), ('JQ_AF_SK', '0.33,0.21')]), tag.format("C,T", "filter", "", format_dict, 0))
#         
        format_dict = OrderedDict(zip("A:TAR:DP2".split(":"), "1:13,32:53".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('TAR', '13,32'), ('DP2', '53'), ('JQ_AF_SK', '0.6')]), tag.format("C", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:TAR:DP2".split(":"), "1:13,50:53".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('TAR', '13,50'), ('DP2', '53'), ('JQ_AF_SK', '0.94')]), tag.format("C", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:TAR:DP2".split(":"), "1:13,70:53".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('TAR', '13,70'), ('DP2', '53'), ('JQ_AF_SK', '1.00')]), tag.format("C", "filter", "", format_dict, 0))

class Strelka_DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_DP_MT,Number=1,Type=Float,Description="Jacquard depth for MuTect (based on DP),Source="Jacquard",Version={0}>\n'.format(__version__), Mutect_DepthTag().metaheader)
                
    def test_format_missingDP2AUTags(self):
        tag = Strelka_DepthTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2')]), tag.format("C", "filter", "", format_dict, 0))
                
    def test_format(self):
        tag = Strelka_DepthTag()
        format_dict = OrderedDict(zip("A:DP2".split(":"), "1:42".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('DP2', '42'), ('JQ_DP_SK', '42')]), tag.format("C", "filter", "", format_dict, 0))
        
        format_dict = OrderedDict(zip("A:AU:CU:GU:TU".split(":"), "1:42,42:5,5:13,13:4,4".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('AU', '42,42'), ('CU', '5,5'), ('GU', '13,13'), ('TU', '4,4'), ('JQ_DP_SK', '64')]), tag.format("C", "filter", "", format_dict, 0))

class Strelka_SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_SOM_MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag),Source="Jacquard",Version={0}>\n'.format(__version__), Mutect_SomaticTag().metaheader)
                
    def test_format_missingPASSTag(self):
        tag = Strelka_SomaticTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        format_dict = OrderedDict(zip(format_param_string.split(":"), format_value_string.split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "reject", "", format_dict, 0))
        self.assertEqual(OrderedDict([('A', '1'), ('B', '2'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "PASS;foo", "", format_dict, 0))
                
    def test_format(self):
        tag = Strelka_SomaticTag()
        format_dict = OrderedDict(zip( "A:SS".split(":"), "1:2".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('SS', '2'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "PASS", "", format_dict, 0))
        self.assertEqual(OrderedDict([('A', '1'), ('SS', '2'), ('JQ_HC_SOM_SK', '1')]), tag.format("C", "PASS", "", format_dict, 1))
        
        format_dict = OrderedDict(zip( "A:SS".split(":"), "1:1".split(":")))
        self.assertEqual(OrderedDict([('A', '1'), ('SS', '1'), ('JQ_HC_SOM_SK', '0')]), tag.format("C", "PASS", "", format_dict, 0))
        self.assertEqual(OrderedDict([('A', '1'), ('SS', '1'), ('JQ_HC_SOM_SK', '1')]), tag.format("C", "PASS", "", format_dict, 1))

class LineProcessorTestCase(unittest.TestCase):
    def test_process_line_singleSample(self):
        tag = MockAddTag("JQ", "42")    
        processor = LineProcessor([tag])
        input_line = "chr1|42|.|ref|alt|qual|filter|INFO|A:B:C|X:Y:Z\n".replace("|", "\t")

        actual_line = processor.add_tags(input_line)
        self.assertEqual(True, actual_line.endswith("\n"))
        self.assertEquals(1, len(actual_line.splitlines()))
        
        trimmed_actual_line = actual_line.rstrip("\n")
        actual_format_params = trimmed_actual_line.split("\t")[8]
        actual_format_values = trimmed_actual_line.split("\t")[9:]        
        self.assertEqual("A:B:C:JQ", actual_format_params)
        self.assertEqual(["X:Y:Z:42"], actual_format_values)
        
    def test_process_line_MultSamples(self):
        tag = MockLowerTag()    
        processor = LineProcessor([tag])
        input_line = "chr1|42|.|ref|alt|qual|filter|INFO|A:B:C|U:V:W|X:Y:Z\n".replace("|", "\t")

        actual_line = processor.add_tags(input_line)
        actual_format_params = actual_line.split("\t")[8]
        actual_format_values = actual_line.split("\t")[9:]
        
        self.assertEqual("a:b:c", actual_format_params)
        self.assertEqual(["u:v:w", "x:y:z\n"], actual_format_values)
        self.assertEqual("chr1\t42\t.\tref\talt\tqual\tfilter\tINFO\ta:b:c\tu:v:w\tx:y:z\n", actual_line)
        
        
    def test_process_line_MultTags(self):
        processor = LineProcessor([MockAddTag("JQ", "42"), MockAddTag("JQ_2", "43"), MockAddTag("JQ_3", "43")])
        input_line = "chr1|42|.|ref|alt|qual|filter|INFO|A:B:C|X:Y:Z\n".replace("|", "\t")

        actual_line = processor.add_tags(input_line)
        self.assertEqual(True, actual_line.endswith("\n"))
        self.assertEquals(1, len(actual_line.splitlines()))
        
        trimmed_actual_line = actual_line.rstrip("\n")
        actual_format_params = trimmed_actual_line.split("\t")[8]
        actual_format_values = trimmed_actual_line.split("\t")[9:]        
        self.assertEqual("A:B:C:JQ:JQ_2:JQ_3", actual_format_params)
        self.assertEqual(["X:Y:Z:42:43:43"], actual_format_values)

    def test_process_line_MultSamples_MultTags(self):
        processor = LineProcessor([MockAddTag("JQ", "42"), MockAddTag("JQ_2", "43"), MockAddTag("JQ_3", "43")])
        input_line = "chr1|42|.|ref|alt|qual|filter|INFO|A:B:C|U:V:W|X:Y:Z\n".replace("|", "\t")

        actual_line = processor.add_tags(input_line)
        self.assertEqual(True, actual_line.endswith("\n"))
        self.assertEquals(1, len(actual_line.splitlines()))
        
        trimmed_actual_line = actual_line.rstrip("\n")
        actual_format_params = trimmed_actual_line.split("\t")[8]
        actual_format_values = trimmed_actual_line.split("\t")[9:]        
        self.assertEqual("A:B:C:JQ:JQ_2:JQ_3", actual_format_params)
        self.assertEqual(["U:V:W:42:43:43","X:Y:Z:42:43:43"], actual_format_values)

class FileProcessorTestCase(unittest.TestCase):
    def test_process_diesIfMissingNormalTumorHeaders(self):
        mockWriter = MockWriter()
        input_metaheaders = ["jacquard.version=X","jacquard.tagMutect.command=foo"]
        processor = FileProcessor(tags=[], execution_context_metadataheaders=input_metaheaders)
        with self.assertRaises(SystemExit) as cm:
            processor.process(reader=["#CHROM\n"], writer=mockWriter, caller="VarScan")
        self.assertEqual(cm.exception.code, 1)
        
    def test_process_prependsExecutionContext(self):
        mockWriter = MockWriter()
        input_metaheaders = ["##jacquard.version=X","##jacquard.tagMutect.command=foo"]
        processor = FileProcessor(tags=[], execution_context_metadataheaders=input_metaheaders)
        processor.process(reader=["#CHROM\tNORMAL\tTUMOR\n"], writer=mockWriter, caller="MuTect")
        actualLines = mockWriter.lines()
        self.assertEqual(3, len(actualLines))
        self.assertEqual("##jacquard.version=X", actualLines[0])
        self.assertEqual("##jacquard.tagMutect.command=foo", actualLines[1])

    def test_prependsExecutionContextWhenBlank(self):
        processor = FileProcessor()
        mockWriter = MockWriter()
        processor.process(reader=[], writer=mockWriter, caller="MuTect")
        self.assertEqual(0, len(mockWriter.lines()))

    def test_process_passthroughExistingMetaHeaders(self):
        mockWriter = MockWriter()
        processor = FileProcessor()
        reader = ["##Hello\n","##World\n"]
        processor.process(reader, mockWriter, "MuTect")
        self.assertEqual(["##Hello", "##World"], mockWriter.lines())

    def test_process_addsNewMetaHeaders(self):
        mockWriter = MockWriter()
        mockTag = MockLowerTag("##mockMetaHeader")
        processor = FileProcessor(tags=[mockTag])
        processor.process(reader=["#CHROM\tNORMAL\tTUMOR\n"], writer=mockWriter, caller="MuTect")
        self.assertEqual(2, len(mockWriter.lines()))
        self.assertEqual("##mockMetaHeader", mockWriter.lines()[0])

    def test_process_adjustsVcfRecords(self):
        mockWriter = MockWriter()
        mockTag = MockLowerTag("##mockMetaHeader")
        recordHeader = "#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR\n".replace("|","\t")
        reader = [ recordHeader,
                  "chr1|1|.|ref|alt|qual|filter|INFO|A:B:C|A1:B1:C1|A2:B2:C2".replace("|","\t"),
                  "chr2|10|.|ref|alt|qual|filter|INFO|A:B:C|A10:B10:C10|A20:B20:C20".replace("|","\t")
                  ]
        processor = FileProcessor(tags=[mockTag])
        processor.process(reader, mockWriter, "MuTect")
        self.assertEqual(4, len(mockWriter.lines()))
        self.assertEqual("##mockMetaHeader", mockWriter.lines()[0])
        self.assertEqual(recordHeader.strip("\n"), mockWriter.lines()[1])
        self.assertEqual("chr1\t1\t.\tref\talt\tqual\tfilter\tINFO\ta:b:c\ta1:b1:c1\ta2:b2:c2", mockWriter.lines()[2])
        self.assertEqual("chr2\t10\t.\tref\talt\tqual\tfilter\tINFO\ta:b:c\ta10:b10:c10\ta20:b20:c20", mockWriter.lines()[3])

class TagVarScanTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        self.saved_stdout = sys.stdout
        sys.stdout = self.output

    def tearDown(self):
        self.output.close()
        sys.stdout = self.saved_stdout
    
    def test_tagVarScanFiles(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("B.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("C.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")

            tag_files(input_dir.path, output_dir.path, [Varscan()])
            
            actual_files = sorted(listdir(output_dir.path))
            
            self.assertEqual(3, len(actual_files))
            self.assertEqual("A_jacquard.vcf", actual_files[0])
            self.assertEqual("B_jacquard.vcf", actual_files[1])
            self.assertEqual("C_jacquard.vcf", actual_files[2])
             
            input_dir.cleanup()
            output_dir.cleanup()

    def test_tagVarScanFiles_content(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("B.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("C.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")

            tag_files(input_dir.path, output_dir.path, [Mutect(), Varscan(), Unknown()])
            actual_files = sorted(listdir(output_dir.path))
            for actual_file in actual_files:
                result = open(output_dir.path + "/" + actual_file).read()
                split_result = result.split("\n")
                self.assertEqual('##source=VarScan2', split_result[0])
                self.assertEqual('##FORMAT=<ID=JQ_AF_VS,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ),Source="Jacquard",Version={0}>'.format(__version__), split_result[1])
                self.assertEqual('##FORMAT=<ID=JQ_DP_VS,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP),Source="Jacquard",Version={0}>'.format(__version__), split_result[2])
                self.assertEqual('##FORMAT=<ID=JQ_SOM_VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1= somatic (based on SOMATIC info tag and if sample is TUMOR),Source="Jacquard",Version={0}>'.format(__version__), split_result[3])
                self.assertEqual('##jacquard.tag.caller=VarScan'.format(__version__), split_result[4])
                self.assertEqual('#CHROM\tNORMAL\tTUMOR', split_result[5])
            
            input_dir.cleanup()
            output_dir.cleanup()

    def test_tagVarScanFiles_ignoresNonVcfFiles(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.txt","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("B.vcf.bak","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("C.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")

            tag_files(input_dir.path, output_dir.path, [Mutect(), Varscan(), Unknown()])
            
            actual_files = sorted(listdir(output_dir.path))
            self.assertEqual(1, len(actual_files))
            
            input_dir.cleanup()
            output_dir.cleanup()
            
    def test_tagVarScanFilestest_tagMutectFiles_ignoresDirectories(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("B.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")

            tag_files(input_dir.path, output_dir.path, [Mutect(), Varscan(), Unknown()], ["execution metaheaders"])
            
            output_list = self.output.getvalue().splitlines()

            self.assertEqual(10, len(output_list))
            self.assertEqual(True, output_list[0].startswith("execution"))
            self.assertEqual(True, output_list[1].startswith("Processing [2]"))
            self.assertEqual(True, output_list[9].startswith("Wrote [2]"))
            
    def test_tagVarScanFilestest_tagMutectFiles_inputDirectoryNoVCFs(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/tag_varscan_test/noVCFs"
        os.mkdir(input_dir)
        output_dir = script_dir + "/tag_varscan_test/output"
        with self.assertRaises(SystemExit) as cm:
            tag_files(input_dir, output_dir, [Mutect(), Varscan(), Unknown()])
        os.rmdir(input_dir)
        self.assertEqual(cm.exception.code, 1)
        
    def test_tagVarScanFiles_noHeaders(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/test_input/unknownCaller"
        output_dir = script_dir + "/test_output"
        with self.assertRaises(SystemExit) as cm:
           tag_files(input_dir, output_dir, [Mutect(), Varscan(), Unknown()])

        self.assertEqual(cm.exception.code, 1)

class  DetermineFileTypesTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        self.saved_stdout = sys.stdout
        sys.stdout = self.output

    def tearDown(self):
        self.output.close()
        sys.stdout = self.saved_stdout
        
    def test_determineFileTypes_withUnknown(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/multi_caller_input/withUnknowns"
        in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
        callers = [Mutect(), Varscan(), Unknown()]
        file_types, inferred_callers = determine_file_types(input_dir, in_files, callers)
        self.assertEqual(["Unknown", "VarScan", "MuTect"], file_types.keys())
        self.assertEqual('tiny_unknown_input.vcf', os.path.basename(file_types.values()[0][0]))
        self.assertEqual('tiny_varscan_input.vcf', os.path.basename(file_types.values()[1][0]))
        self.assertEqual('tiny_mutect_input.vcf', os.path.basename(file_types.values()[2][0]))
        self.assertEqual('tiny_mutect_input2.vcf', os.path.basename(file_types.values()[2][1]))
        
        with self.assertRaises(SystemExit) as cm:
            print_file_types(file_types)
        self.assertEqual(cm.exception.code, 1)
        
        output_list = self.output.getvalue().splitlines()
                
        self.assertEqual("tiny_mutect_input.vcf: ##jacquard.tag.handler=MuTect", output_list[0])
        self.assertEqual("tiny_mutect_input2.vcf: ##jacquard.tag.handler=MuTect", output_list[2])
        self.assertEqual("ERROR: tiny_unknown_input.vcf: ##jacquard.tag.handler=Unknown", output_list[4])
        self.assertEqual("tiny_varscan_input.vcf: ##jacquard.tag.handler=VarScan", output_list[5])
        self.assertEqual("Recognized [1] Unknown file(s)", output_list[7])
        self.assertEqual("Recognized [1] VarScan file(s)", output_list[8])
        self.assertEqual("Recognized [2] MuTect file(s)", output_list[9])
        
        
    def test_determineFileTypes_noUnknown(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/multi_caller_input/"
        in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
        
        callers = [Mutect(), Varscan(), Unknown()]
        file_types, inferred_callers = determine_file_types(input_dir, in_files, callers)
        self.assertEqual(["VarScan", "MuTect"], file_types.keys())
        self.assertEqual('tiny_varscan_input.vcf', os.path.basename(file_types.values()[0][0]))
        self.assertEqual('tiny_mutect_input.vcf', os.path.basename(file_types.values()[1][0]))
        self.assertEqual('tiny_mutect_input2.vcf', os.path.basename(file_types.values()[1][1]))

        print_file_types(file_types)
        output_list = self.output.getvalue().splitlines()
        
        self.assertEqual("tiny_mutect_input.vcf: ##jacquard.tag.handler=MuTect", output_list[0])
        self.assertEqual("tiny_mutect_input2.vcf: ##jacquard.tag.handler=MuTect", output_list[2])
        self.assertEqual("tiny_varscan_input.vcf: ##jacquard.tag.handler=VarScan", output_list[4])
        self.assertEqual("Recognized [1] VarScan file(s)", output_list[6])
        self.assertEqual("Recognized [2] MuTect file(s)", output_list[7])


class MockLowerTag():
    def __init__(self, metaheader=""):
        self.metaheader = metaheader

    def format(self, alt, filter, info, format_dict, count):
        str_keys = ",".join(format_dict.keys())
        new_keys = str_keys.lower()
        
        str_vals = ",".join(format_dict.values())
        new_vals = str_vals.lower()
        
        new_dict = OrderedDict(zip(new_keys.split(","), new_vals.split(",")))
        
        return new_dict
    
class MockAddTag():
    def __init__(self, tag="MockFormat", value="MockValue"):
        self.metaheader = "MockMetaheader"
        self.tag= tag
        self.value = value

    def format(self, alt, filter, info, format_dict, count):
        new_keys = format_dict.keys()
        new_keys.append(self.tag)
        new_vals = format_dict.values()
        new_vals.append(self.value)

        new_dict = OrderedDict(zip(new_keys, new_vals))
        
        return new_dict

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
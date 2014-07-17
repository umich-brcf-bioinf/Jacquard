#!/usr/bin/python2.7
import os
from os import listdir
from StringIO import StringIO
import sys
from stat import *
import testfixtures
from testfixtures import TempDirectory
import unittest
from bin.tag_varscan import AlleleFreqTag, DepthTag, SomaticTag, LineProcessor, FileProcessor, tag_varscan_files, validate_directories

class AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_VS,Number=1,Type=Float, Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ).">\n', AlleleFreqTag().metaheader)
                
    def test_format_missingAFTag(self):
        tag = AlleleFreqTag()
        info_string = ""
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(info_string, format_param_string, format_value_string))
                
    def test_format_rounds(self):
        tag = AlleleFreqTag()
        self.assertEqual(("A:FREQ:JQ_AF_VS", "1:0.2:0.2"), tag.format("", "A:FREQ", "1:0.2"))
        self.assertEqual(("A:FREQ:JQ_AF_VS", "1:0.20:0.20"), tag.format("", "A:FREQ", "1:0.20"))
        self.assertEqual(("A:FREQ:JQ_AF_VS", "1:0.204:0.2"), tag.format("", "A:FREQ", "1:0.204"))
        self.assertEqual(("A:FREQ:JQ_AF_VS", "1:0.205:0.21"), tag.format("", "A:FREQ", "1:0.205"))
        self.assertEqual(("A:FREQ:JQ_AF_VS", "1:0.206:0.21"), tag.format("", "A:FREQ", "1:0.206"))
        self.assertEqual(("A:FREQ:JQ_AF_VS", "1:1.0:1.0"), tag.format("", "A:FREQ", "1:1.0"))
        self.assertEqual(("A:FREQ:JQ_AF_VS", "1:1.00:1.00"), tag.format("", "A:FREQ", "1:1.00"))

class DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_DP_VS,Number=1,Type=Float, Description="Jacquard depth for VarScan (based on DP).">\n', DepthTag().metaheader)
                
    def test_format_missingDPTag(self):
        tag = DepthTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format("", format_param_string, format_value_string))
                
    def test_format(self):
        tag = DepthTag()
        self.assertEqual(("A:DP:JQ_DP_VS", "1:42:42"), tag.format("", "A:DP", "1:42"))

class SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_SOM_VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1= somatic (based on SOMATIC info tag and if sample is TUMOR).">\n', SomaticTag().metaheader)
                
    def test_format_missingSSInfoTag(self):
        tag = SomaticTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format("INFO", format_param_string, format_value_string))
                
    def test_format(self):
        tag = SomaticTag()
        self.assertEqual(("A:JQ_SOM_VS", "1:1"), tag.format("INFO;SS=2", "A", "1"))
        self.assertEqual(("A", "1"), tag.format("INFO", "A", "1"))

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
    
    def test_process_prependsExecutionContext(self):
        mockWriter = MockWriter()
        input_metaheaders = ["jacquard.version=X","jacquard.tagMutect.command=foo"]
        processor = FileProcessor(tags=[], execution_context_metadataheaders=input_metaheaders)
        processor.process(reader=["#CHROM\n"], writer=mockWriter)
        actualLines = mockWriter.lines()
        self.assertEqual(3, len(actualLines))
        self.assertEqual("##jacquard.version=X", actualLines[0])
        self.assertEqual("##jacquard.tagMutect.command=foo", actualLines[1])

    def test_prependsExecutionContextWhenBlank(self):
        processor = FileProcessor()
        mockWriter = MockWriter()
        processor.process(reader=[], writer=mockWriter)
        self.assertEqual(0, len(mockWriter.lines()))

    def test_process_passthroughExistingMetaHeaders(self):
        mockWriter = MockWriter()
        processor = FileProcessor()
        reader = ["##Hello\n","##World\n"]
        processor.process(reader, mockWriter)
        self.assertEqual(["##Hello", "##World"], mockWriter.lines())

    def test_process_addsNewMetaHeaders(self):
        mockWriter = MockWriter()
        mockTag = MockLowerTag("##mockMetaHeader")
        processor = FileProcessor(tags=[mockTag])
        processor.process(reader=["#CHROM\n"], writer=mockWriter)
        self.assertEqual(2, len(mockWriter.lines()))
        self.assertEqual("##mockMetaHeader", mockWriter.lines()[0])

    def test_process_adjustsVcfRecords(self):
        mockWriter = MockWriter()
        mockTag = MockLowerTag("##mockMetaHeader")
        recordHeader = "#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleA|SampleB".replace("|","\t")
        reader = [ recordHeader,
                  "chr1|1|.|ref|alt|qual|filter|INFO|A:B:C|A1:B1:C1|A2:B2:C2".replace("|","\t"),
                  "chr2|10|.|ref|alt|qual|filter|INFO|A:B:C|A10:B10:C10|A20:B20:C20".replace("|","\t")
                  ]
        processor = FileProcessor(tags=[mockTag])
        processor.process(reader, mockWriter)
        self.assertEqual(4, len(mockWriter.lines()))
        self.assertEqual("##mockMetaHeader", mockWriter.lines()[0])
        self.assertEqual(recordHeader, mockWriter.lines()[1])
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
            input_dir.write("A.vcf","##source=VarScan2\n#CHROM\n")
            input_dir.write("B.vcf","##source=VarScan2\n#CHROM\n")
            input_dir.write("C.vcf","##source=VarScan2\n#CHROM\n")

            tag_varscan_files(input_dir.path, output_dir.path)
            actual_files = sorted(listdir(output_dir.path))

            self.assertEqual(3, len(actual_files))
            self.assertEqual("A_jacquard.vcf", actual_files[0])
            self.assertEqual("B_jacquard.vcf", actual_files[1])
            self.assertEqual("C_jacquard.vcf", actual_files[2])
             
            input_dir.cleanup()
            output_dir.cleanup()

    def test_tagVarScanFiles_content(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","##source=VarScan2\n#CHROM\n")
            input_dir.write("B.vcf","##source=VarScan2\n#CHROM\n")
            input_dir.write("C.vcf","##source=VarScan2\n#CHROM\n")

            tag_varscan_files(input_dir.path, output_dir.path)
            actual_files = sorted(listdir(output_dir.path))
            for actual_file in actual_files:
                result = open(output_dir.path + "/" + actual_file).read()
                split_result = result.split("\n")
                self.assertEqual('##source=VarScan2', split_result[0])
                self.assertEqual('##FORMAT=<ID=JQ_AF_VS,Number=1,Type=Float, Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ).">', split_result[1])
                self.assertEqual('##FORMAT=<ID=JQ_DP_VS,Number=1,Type=Float, Description="Jacquard depth for VarScan (based on DP).">', split_result[2])
                self.assertEqual('##FORMAT=<ID=JQ_SOM_VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1= somatic (based on SOMATIC info tag and if sample is TUMOR).">', split_result[3])
                self.assertEqual('#CHROM', split_result[4])
            
            input_dir.cleanup()
            output_dir.cleanup()

    def test_tagVarScanFiles_noMutectHeader(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","#CHROM\n")
            input_dir.write("B.vcf","#CHROM\n")
            input_dir.write("C.vcf","#CHROM\n")

            tag_varscan_files(input_dir.path, output_dir.path)
            actual_files = sorted(listdir(output_dir.path))
            input_files = sorted(listdir(input_dir.path))
            
            self.assertEqual(0, len(actual_files))
            self.assertEqual(3, len(input_files))
             
            input_dir.cleanup()
            output_dir.cleanup()

    def test_tagVarScanFiles_ignoresNonVcfFiles(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.txt","##source=VarScan2\n#CHROM\n")
            input_dir.write("B.vcf.bak","##source=VarScan2\n#CHROM\n")
            input_dir.write("C.vcf","##source=VarScan2\n#CHROM\n")

            tag_varscan_files(input_dir.path, output_dir.path)
            
            actual_files = sorted(listdir(output_dir.path))
            self.assertEqual(1, len(actual_files))
            
            input_dir.cleanup()
            output_dir.cleanup()

    def test_tagVarScanFiles_ignoresDirectories(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.makedir("A")
            input_dir.makedir("B.vcf")

            tag_varscan_files(input_dir.path, output_dir.path)
            
            actual_files = sorted(listdir(output_dir.path))
            self.assertEqual(0, len(actual_files))
            
            input_dir.cleanup()
            output_dir.cleanup()
            
    def test_tagVarScanFilestest_tagMutectFiles_ignoresDirectories(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","##source=VarScan2\n#CHROM\n")
            input_dir.write("B.vcf","#CHROM\n")

            tag_varscan_files(input_dir.path, output_dir.path, ["execution metaheaders"])
            
            output_list = self.output.getvalue().splitlines()
            self.assertEqual(3, len(output_list))
            self.assertEqual(True, output_list[0].startswith("execution"))
            self.assertEqual(True, output_list[1].startswith("Processing [2]"))
            self.assertEqual(True, output_list[2].startswith("Wrote [1]"))
            
    def test_tagVarScanFilestest_tagMutectFiles_inputDirectoryNoVCFs(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/tag_varscan_test/noVCFs"
        os.mkdir(input_dir)
        output_dir = script_dir + "/tag_varscan_test/output"
        with self.assertRaises(SystemExit) as cm:
            tag_varscan_files(input_dir, output_dir)
        os.rmdir(input_dir)
        self.assertEqual(cm.exception.code, 1)
        

class test_ValidateDirectoriesTestCase(unittest.TestCase):
    def test_validateDirectories_inputDirectoryDoesntExist(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/tag_varscan_test/foo"
        output_dir = script_dir + "/tag_varscan_test/output"
        
        with self.assertRaises(SystemExit) as cm:
            validate_directories(input_dir, output_dir)
        self.assertEqual(cm.exception.code, 1)
    
    def test_validateDirectories_inputDirectoryUnreadable(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/tag_varscan_test/unreadable"
        output_dir = script_dir + "/tag_varscan_test/output"

        with self.assertRaises(SystemExit) as cm:
            validate_directories(input_dir, output_dir)
        self.assertEqual(cm.exception.code, 1)
        
    def test_validateDirectories_outputDirectoryNotCreated(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/tag_varscan_test/input"
        first_out_dir = script_dir + "/tag_varscan_test/unwriteable"
        
        with self.assertRaises(SystemExit) as cm:
            validate_directories(input_dir, first_out_dir + "/foo")
        self.assertEqual(cm.exception.code, 1)

class MockLowerTag():
    def __init__(self, metaheader=""):
        self.metaheader = metaheader

    def format(self, info, params, values):
        return (params.lower(), values.lower())
    
class MockAddTag():
    def __init__(self, tag="MockFormat", value="MockValue"):
        self.metaheader = "MockMetaheader"
        self.tag= tag
        self.value = value

    def format(self, info, params, values):
        return (params+":"+self.tag, values+":"+self.value)

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
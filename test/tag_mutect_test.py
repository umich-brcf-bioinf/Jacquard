#!/usr/bin/python2.7
import os
from os import listdir
import testfixtures
from testfixtures import TempDirectory
import unittest
from bin.tag_mutect import AlleleFreqTag, DepthTag, SomaticTag, LineProcessor, FileProcessor, tag_mutect_files

class AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_MT,Number=1,Type=Float, Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA).">\n', AlleleFreqTag().metaheader)
                
    def test_format_missingAFTag(self):
        tag = AlleleFreqTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(format_param_string, format_value_string))
                
    def test_format_rounds(self):
        tag = AlleleFreqTag()
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.2:0.2"), tag.format("A:FA", "1:0.2"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.20:0.20"), tag.format("A:FA", "1:0.20"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.204:0.2"), tag.format("A:FA", "1:0.204"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.205:0.21"), tag.format("A:FA", "1:0.205"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.206:0.21"), tag.format("A:FA", "1:0.206"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:1.0:1.0"), tag.format("A:FA", "1:1.0"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:1.00:1.00"), tag.format("A:FA", "1:1.00"))

class DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_DP_MT,Number=1,Type=Float, Description="Jacquard depth for MuTect (based on DP).">\n', DepthTag().metaheader)
                
    def test_format_missingDPTag(self):
        tag = DepthTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(format_param_string, format_value_string))
                
    def test_format(self):
        tag = DepthTag()
        self.assertEqual(("A:DP:JQ_DP_MT", "1:42:42"), tag.format("A:DP", "1:42"))

class SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_SOM_MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag).">\n', SomaticTag().metaheader)
                
    def test_format_missingSSTag(self):
        tag = SomaticTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(format_param_string, format_value_string))
                
    def test_format(self):
        tag = SomaticTag()
        self.assertEqual(("A:SS:JQ_SOM_MT", "1:2:1"), tag.format("A:SS", "1:2"))
        self.assertEqual(("A:SS:JQ_SOM_MT", "1:1:0"), tag.format("A:SS", "1:1"))

class LineProcessorTestCase(unittest.TestCase):
    def test_process_line_singleSample(self):
        tag = MockLowerTag()    
        processor = LineProcessor([tag])
        input_line = "chr1|42|.|ref|alt|qual|filter|INFO|A:B:C|X:Y:Z".replace("|", "\t")

        actual_line = processor.add_tags(input_line)
        actual_format_params = actual_line.split("\t")[8]
        actual_format_values = actual_line.split("\t")[9:]
        
        self.assertEqual("a:b:c", actual_format_params)
        self.assertEqual(["x:y:z"], actual_format_values)
        self.assertEqual("chr1\t42\t.\tref\talt\tqual\tfilter\tINFO\ta:b:c\tx:y:z", actual_line)
        
    def test_process_line_MultSample(self):
        tag = MockLowerTag()    
        processor = LineProcessor([tag])
        input_line = "chr1|42|.|ref|alt|qual|filter|INFO|A:B:C|U:V:W|X:Y:Z".replace("|", "\t")

        actual_line = processor.add_tags(input_line)
        actual_format_params = actual_line.split("\t")[8]
        actual_format_values = actual_line.split("\t")[9:]
        
        self.assertEqual("a:b:c", actual_format_params)
        self.assertEqual(["u:v:w", "x:y:z"], actual_format_values)
        self.assertEqual("chr1\t42\t.\tref\talt\tqual\tfilter\tINFO\ta:b:c\tu:v:w\tx:y:z", actual_line)

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


class TagMutectTest(unittest.TestCase):
    def test_tagMutectFiles(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","##MuTect\n#CHROM\n")
            input_dir.write("B.vcf","##MuTect\n#CHROM\n")
            input_dir.write("C.vcf","##MuTect\n#CHROM\n")

            tag_mutect_files(input_dir.path, output_dir.path)
            actual_files = sorted(listdir(output_dir.path))

            self.assertEqual(3, len(actual_files))
            self.assertEqual("A_jacquard.vcf", actual_files[0])
            self.assertEqual("B_jacquard.vcf", actual_files[1])
            self.assertEqual("C_jacquard.vcf", actual_files[2])
             
            input_dir.cleanup()
            output_dir.cleanup()

    def test_tagMutectFilesContent(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","##MuTect\n#CHROM\n")
            input_dir.write("B.vcf","##MuTect\n#CHROM\n")
            input_dir.write("C.vcf","##MuTect\n#CHROM\n")

            tag_mutect_files(input_dir.path, output_dir.path)
            actual_files = sorted(listdir(output_dir.path))
            for actual_file in actual_files:
                result = open(output_dir.path + "/" + actual_file).read()
                self.assertEqual('##MuTect\n##FORMAT=<ID=JQ_AF_MT,Number=1,Type=Float, Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA).">\n#CHROM\n', result)
            
            input_dir.cleanup()
            output_dir.cleanup()

    def test_tagMutectFilesNoMutectHeader(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","#CHROM\n")
            input_dir.write("B.vcf","#CHROM\n")
            input_dir.write("C.vcf","#CHROM\n")

            tag_mutect_files(input_dir.path, output_dir.path)
            actual_files = sorted(listdir(output_dir.path))
            input_files = sorted(listdir(input_dir.path))
            
            self.assertEqual(0, len(actual_files))
            self.assertEqual(3, len(input_files))
             
            input_dir.cleanup()
            output_dir.cleanup()
    
#     def test_tagMutectFiles(self):
#         #make a mutect file
#         script_dir = os.path.dirname(os.path.abspath(__file__))
#         input_path = script_dir + "/test_mutect/input"
#         expected_output = script_dir + "/test_mutect/expected/expected_tiny_mutect_output.vcf"
#        
#         with TempDirectory() as d:
#             process_mutect_file(input_path, d)
#             actual_output = output_path + "/jacquard_tiny_mutect_output.vcf"
#             
#             self.compare(d.read(expected_output), actual_output)
        



class MockLowerTag():
    def __init__(self, metaheader=""):
        self.metaheader = metaheader

    def format(self, params, values):
        return (params.lower(), values.lower())
    
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
from collections import OrderedDict
import glob
import os
from os import listdir
from StringIO import StringIO
import sys
from testfixtures import TempDirectory
import unittest
from jacquard.tag import LineProcessor, FileProcessor, tag_files, determine_file_types, print_file_types
from jacquard.variant_callers import varscan, mutect, unknown
from jacquard.jacquard_utils import __version__


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
        script_dir = os.path.dirname(os.path.abspath(__file__))
        output_dir = script_dir+ "/reference_files/fake_output/"
        try:
            os.mkdir(output_dir)
        except:
            pass
        
        input_metaheaders = ["jacquard.version=X","jacquard.tagMutect.command=foo"]
        processor = FileProcessor(output_dir, tags=[], execution_context_metadataheaders=input_metaheaders)
        with self.assertRaises(SystemExit) as cm:
            processor.process(reader=["#CHROM\n"], writer=mockWriter, in_file_name="foo", caller="VarScan", handlers={"foo":""})
        self.assertEqual(cm.exception.code, 1)
        
    def test_process_prependsExecutionContext(self):
        mockWriter = MockWriter()
        input_metaheaders = ["##jacquard.version=X","##jacquard.tagMutect.command=foo"]
        processor = FileProcessor("foo_output_dir", tags=[], execution_context_metadataheaders=input_metaheaders)
        processor.process(reader=["#CHROM\tNORMAL\tTUMOR\n"], writer=mockWriter, in_file_name="foo", caller="MuTect", handlers={"foo":""})
        actualLines = mockWriter.lines()
        self.assertEqual(5, len(actualLines))
        self.assertEqual("##jacquard.version=X", actualLines[1])
        self.assertEqual("##jacquard.tagMutect.command=foo", actualLines[2])

    def test_prependsExecutionContextWhenBlank(self):
        processor = FileProcessor("foo_output_dir")
        mockWriter = MockWriter()
        processor.process(reader=[], writer=mockWriter, in_file_name="foo", caller="MuTect", handlers={"foo":""})
        self.assertEqual(1, len(mockWriter.lines()))

    def test_process_passthroughExistingMetaHeaders(self):
        mockWriter = MockWriter()
        processor = FileProcessor("foo_output_dir")
        reader = ["##Hello\n","##World\n"]
        processor.process(reader, mockWriter, "foo", "MuTect", handlers={"foo":""})
        self.assertEqual(["", "##Hello", "##World"], mockWriter.lines())

    def test_process_addsNewMetaHeaders(self):
        mockWriter = MockWriter()
        mockTag = MockLowerTag("##mockMetaHeader")
        processor = FileProcessor("foo_output_dir", tags=[mockTag])
        processor.process(reader=["#CHROM\tNORMAL\tTUMOR\n"], writer=mockWriter, in_file_name="foo", caller="MuTect", handlers={"foo":""})
        self.assertEqual(4, len(mockWriter.lines()))
        self.assertEqual("##mockMetaHeader", mockWriter.lines()[1])

    def test_process_adjustsVcfRecords(self):
        mockWriter = MockWriter()
        mockTag = MockLowerTag("##mockMetaHeader")
        recordHeader = "#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR\n".replace("|","\t")
        reader = [ recordHeader,
                  "chr1|1|.|ref|alt|qual|filter|INFO|A:B:C|A1:B1:C1|A2:B2:C2".replace("|","\t"),
                  "chr2|10|.|ref|alt|qual|filter|INFO|A:B:C|A10:B10:C10|A20:B20:C20".replace("|","\t")
                  ]
        processor = FileProcessor("foo_output_dir", tags=[mockTag])
        processor.process(reader, mockWriter, "foo", "MuTect", {"foo":""})
        
        self.assertEqual(6, len(mockWriter.lines()))
        self.assertEqual("##mockMetaHeader", mockWriter.lines()[1])
        self.assertEqual("##jacquard.tag.caller=MuTect", mockWriter.lines()[2])
        self.assertEqual(recordHeader.strip("\n"), mockWriter.lines()[3])
        self.assertEqual("chr1\t1\t.\tref\talt\tqual\tfilter\tINFO\ta:b:c\ta1:b1:c1\ta2:b2:c2", mockWriter.lines()[4])
        self.assertEqual("chr2\t10\t.\tref\talt\tqual\tfilter\tINFO\ta:b:c\ta10:b10:c10\ta20:b20:c20", mockWriter.lines()[5])

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


            tag_files(input_dir.path, output_dir.path, [varscan.Varscan()])
            
            actual_files = sorted(listdir(output_dir.path))
            
            self.assertEqual(3, len(actual_files))
            self.assertEqual("A.jacquardTags.vcf", actual_files[0])
            self.assertEqual("B.jacquardTags.vcf", actual_files[1])
            self.assertEqual("C.jacquardTags.vcf", actual_files[2])
             
            input_dir.cleanup()
            output_dir.cleanup()

    def test_tagVarScanFiles_content(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("B.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("C.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")


            tag_files(input_dir.path, output_dir.path, [mutect.Mutect(), varscan.Varscan(), unknown.Unknown()])
            actual_files = sorted(listdir(output_dir.path))
            for actual_file in actual_files:
                result = open(output_dir.path + "/" + actual_file).read()
                split_result = result.split("\n")
                self.assertEqual('##jacquard.tag.handler=VarScan', split_result[0])
                self.assertEqual('##source=VarScan2', split_result[1])
                self.assertEqual('##FORMAT=<ID=JQ_AF_VS,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ),Source="Jacquard",Version={0}>'.format(__version__), split_result[2])
                self.assertEqual('##FORMAT=<ID=JQ_DP_VS,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP),Source="Jacquard",Version={0}>'.format(__version__), split_result[3])
                self.assertEqual('##FORMAT=<ID=JQ_HC_SOM_VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1= somatic (based on SOMATIC info tag and if sample is TUMOR),Source="Jacquard",Version={0}>'.format(__version__), split_result[4])
                self.assertEqual('##jacquard.tag.caller=VarScan'.format(__version__), split_result[5])
                self.assertEqual('#CHROM\tNORMAL\tTUMOR', split_result[6])
            
            input_dir.cleanup()
            output_dir.cleanup()

    def test_tagVarScanFiles_ignoresNonVcfFiles(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.txt","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("B.vcf.bak","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("C.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")


            tag_files(input_dir.path, output_dir.path, [mutect.Mutect(), varscan.Varscan(), unknown.Unknown()])
            
            actual_files = sorted(listdir(output_dir.path))
            self.assertEqual(1, len(actual_files))
            
            input_dir.cleanup()
            output_dir.cleanup()
            
    def test_tagVarScanFilestest_tagMutectFiles_ignoresDirectories(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            input_dir.write("B.vcf","##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")

            tag_files(input_dir.path, output_dir.path, [mutect.Mutect(), varscan.Varscan(), unknown.Unknown()], ["execution metaheaders"])
            output_list = self.output.getvalue().splitlines()
 
            self.assertEqual(10, len(output_list))
            self.assertEqual(True, output_list[0].startswith("execution"))
            self.assertEqual(True, output_list[1].startswith("Processing [2]"))
            self.assertEqual(True, output_list[9].startswith("Wrote [2]"))
            input_dir.cleanup()
            output_dir.cleanup()
            
    def test_tagVarScanFilestest_tagMutectFiles_inputDirectoryNoVCFs(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/reference_files/tag_varscan_test/noVCFs"

        output_dir = script_dir + "/reference_files/tag_varscan_test/output"
        os.mkdir(output_dir)
        with self.assertRaises(SystemExit) as cm:

            tag_files(input_dir, output_dir, [mutect.Mutect(), varscan.Varscan(), unknown.Unknown()])
        
        self.assertEqual(cm.exception.code, 1)
        
    def test_tagVarScanFiles_noHeaders(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/reference_files/test_input/unknownCaller"
        output_dir = script_dir + "/reference_files/invalid_output"
        os.mkdir(output_dir)
        with self.assertRaises(SystemExit) as cm:

            tag_files(input_dir, output_dir, [mutect.Mutect(), varscan.Varscan(), unknown.Unknown()])

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
        input_dir = script_dir + "/reference_files/multi_caller_input/withUnknowns"
        output_dir = script_dir + "/reference_files/multi_caller_output/"
        try:
            os.mkdir(output_dir)
        except:
            pass
        
        in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))

        callers = [mutect.Mutect(), varscan.Varscan(), unknown.Unknown()]
        file_types, handlers, inferred_callers = determine_file_types(input_dir, in_files, callers)
        self.assertEqual(["Unknown", "VarScan", "MuTect"], file_types.keys())
        self.assertEqual('tiny_unknown_input.vcf', os.path.basename(file_types.values()[0][0]))
        self.assertEqual('tiny_varscan_input.vcf', os.path.basename(file_types.values()[1][0]))
        self.assertEqual('tiny_mutect_input.vcf', os.path.basename(file_types.values()[2][0]))
        self.assertEqual('tiny_mutect_input2.vcf', os.path.basename(file_types.values()[2][1]))
        
        with self.assertRaises(SystemExit) as cm:
            print_file_types(output_dir, file_types)
        self.assertEqual(cm.exception.code, 1)
        
        output_list = self.output.getvalue().splitlines()
                
        self.assertEqual("tiny_mutect_input.vcf: ##jacquard.tag.handler=MuTect", output_list[0])
        self.assertEqual("tiny_mutect_input2.vcf: ##jacquard.tag.handler=MuTect", output_list[2])
        self.assertEqual("ERROR. tiny_unknown_input.vcf: ##jacquard.tag.handler=Unknown", output_list[4])
        self.assertEqual("tiny_varscan_input.vcf: ##jacquard.tag.handler=VarScan", output_list[5])
        self.assertEqual("Recognized [1] Unknown file(s)", output_list[7])
        self.assertEqual("Recognized [1] VarScan file(s)", output_list[8])
        self.assertEqual("Recognized [2] MuTect file(s)", output_list[9])
        
        
    def test_determineFileTypes_noUnknown(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/reference_files/multi_caller_input/"
        output_dir = script_dir + "/reference_files/multi_caller_output/"
        in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
        

        callers = [mutect.Mutect(), varscan.Varscan(), unknown.Unknown()]
        file_types, handlers, inferred_callers = determine_file_types(input_dir, in_files, callers)
        self.assertEqual(["VarScan", "MuTect"], file_types.keys())
        self.assertEqual('tiny_varscan_input.vcf', os.path.basename(file_types.values()[0][0]))
        self.assertEqual('tiny_mutect_input.vcf', os.path.basename(file_types.values()[1][0]))
        self.assertEqual('tiny_mutect_input2.vcf', os.path.basename(file_types.values()[1][1]))

        print_file_types(output_dir, file_types)
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
# pylint: disable=C0103,C0301,R0903,R0904
from StringIO import StringIO
from argparse import Namespace

from jacquard.variant_callers import varscan, strelka, variant_caller_factory
from testfixtures import TempDirectory

from jacquard.normalize import _partition_input_files, _determine_caller_per_directory

from jacquard.vcf import FileReader, FileWriter
import jacquard.utils as utils
    
import jacquard.normalize as normalize
import jacquard.utils as utils
import os
import sys
import unittest

# def build_mock_get_caller_method(caller):
#     def get_caller(metaheaders, column_header, name):
#         return caller
#     return get_caller

class MockCallerFactory(object):
    def __init__(self, caller):
        self.caller = caller
        self.last_filename = None

    def get_caller(self, metaheaders, column_header, name):
        self.last_filename = name
        return self.caller
    
class MockCaller(object):
    def __init__(self, name="MockCaller", metaheaders=["##mockMetaheader1"]):
        self.name = name
        self.metaheaders = metaheaders
        self.file_name_search = "snps|indels"

    def add_tags(self, vcfRecord):
        return vcfRecord

    def decorate_files(self, filenames, decorator):
        return filenames[0]+"foo"

    def get_new_metaheaders(self):
        return self.metaheaders

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
    
class MockFileWriter():
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

class NormalizeTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr
        unittest.TestCase.tearDown(self)

    def test_execute(self):
        vcf_content1 = ('''##source=strelka
##file1
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
''').replace('|', "\t")
        vcf_content2 = ('''##source=strelka
##file2
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|10|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|10|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
''').replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("P1.strelka.snvs.vcf", vcf_content1)
            input_dir.write("P1.strelka.indels.vcf", vcf_content2)
            args = Namespace(input_dir=input_dir.path,
                             output_dir=output_dir.path)

            normalize.execute(args, ["extra_header1", "extra_header2"])

            output_dir.check("P1.strelka.normalized.vcf")
            with open(os.path.join(output_dir.path, "P1.strelka.normalized.vcf")) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()

        self.assertEquals(8, len(actual_output_lines), "normalize output wrong number of lines")

    def test_validate_single_caller(self):
        with TempDirectory() as input_dir:
            strelka_file1 = input_dir.write("strelka1.vcf",
                                            "##source=strelka\n#colHeader")
            strelka_file2 = input_dir.write("strelka2.vcf",
                                            "##source=strelka\n#colHeader")
            caller = normalize._validate_single_caller([strelka_file1, strelka_file2],
                                                       variant_caller_factory.get_caller)

            self.assertEquals("Strelka", caller.name)

    def test_validate_single_caller_raisesWhenDistinctCallers(self):
        with TempDirectory() as input_dir:
            strelka_file = input_dir.write("strelka1.vcf",
                                           "##source=strelka\n#colHeader")
            varscan_file = input_dir.write("varscan.vcf",
                                           "##source=VarScan2\n" +\
                                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR")
            self.assertRaises(utils.JQException,
                              normalize._validate_single_caller,
                              [strelka_file, varscan_file],
                              variant_caller_factory.get_caller)

    def test_validate_single_caller_raisesWhenUnrecognizedCaller(self):
        with TempDirectory() as input_dir:
            strelka_file = input_dir.write("strelka1.vcf",
                                           "##source=strelka\n#colHeader")
            unrecognized_file = input_dir.write("foo.vcf",
                                                "##source=foo\n#colHeader")
            self.assertRaises(utils.JQException,
                              normalize._validate_single_caller,
                              [strelka_file, unrecognized_file],
                              variant_caller_factory.get_caller)

    def Ytest__partition_input_files(self):
        in_files = ["A.1.snps.vcf", "A.1.indels.vcf", "B.snps.vcf"]
        output_dir_path = "output_dir_path"
        caller = MockCaller()
        writer_to_readers = _partition_input_files(in_files, output_dir_path, caller)
        
        writerA = FileWriter(os.path.join(output_dir_path,"A.1.normalized.vcf"))
        readersA = [FileReader(os.path.join("A.1.snps.vcf")), 
                    FileReader(os.path.join("A.1.indels.vcf"))]
        writerB = FileWriter(os.path.join(output_dir_path,"B.normalized.vcf"))           
        readersB = [FileReader(os.path.join("B.snps.vcf"))]
        self.assertEquals({writerA: readersA, writerB: readersB}, 
                          writer_to_readers)
        
    def test__partition_input_files(self):
        in_files = ["A.","A.","B."]
        caller = MockCaller()
        output_dir_path = ""
        writer_to_readers = _partition_input_files(in_files, output_dir_path, caller)
        self.maxDiff=None
        writerA = FileWriter("A.foo")
        readersA = [FileReader("A."), 
                    FileReader("A.")]
        writerB = FileWriter("B.foo")           
        readersB = [FileReader("B.")]
        
        self.assertEquals({writerA: readersA, writerB: readersB}, 
                          writer_to_readers)
        
    def test__determine_caller_per_directory(self):
            with TempDirectory() as input_dir:
                A = input_dir.write("A.vcf","##source=strelka\n#colHeader")
                B = input_dir.write("B.vcf","##source=strelka\n#colHeader")
                input_files = [A, B]
                
                mock_caller = MockCaller()
                mock_caller_factory = MockCallerFactory(mock_caller)
                
                caller = _determine_caller_per_directory(input_files, mock_caller_factory.get_caller)
                
                self.assertEquals(mock_caller,caller)
                self.assertEquals(mock_caller_factory.last_filename, "A.vcf")
        
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

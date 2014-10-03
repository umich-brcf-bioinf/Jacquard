# pylint: disable=C0103,C0301,R0903,R0904
from collections import OrderedDict
import glob
import os
from os import listdir
from StringIO import StringIO
import sys
from testfixtures import TempDirectory
import unittest
from jacquard.tag import tag_files
from jacquard.variant_callers import varscan, mutect
from jacquard.jacquard_utils import __version__
import jacquard.tag as tag
import jacquard.jacquard_utils as jacquard_utils
from argparse import Namespace


class MockCallerFactory(object):
    def __init__(self, vcf_to_caller):
        self.vcf_to_caller = vcf_to_caller
        self.last_vcf = None

    def get_vcf_caller(self, vcf):
        self.last_vcf = vcf
        for known_vcf, caller in self.vcf_to_caller.iteritems():
            if known_vcf == vcf:
                return caller
        raise jacquard_utils.JQException("No caller found")
        
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
        
class MockCaller(object):
    def __init__(self, name="MockCaller", metaheaders=["##mockMetaheader1"]):
        self.name = name
        self.metaheaders = metaheaders

    def add_tags(self, vcfRecord):
        return vcfRecord

    def get_new_metaheaders(self):
        return self.metaheaders

class MockVcfReader(object):
    def __init__(self, name="vcfName", metaheaders=["##metaheaders"], column_header="#header"):
        self.name = name
        self.metaheaders = metaheaders
        self.caller = MockCaller()
        self.column_header = column_header
        self.opened = False
        self.closed = False

    def open(self):
        self.opened = True

    def vcf_records(self):
        return iter(["foo"])
    
    def close(self):
        self.closed = True

def build_mock_get_caller_method(caller):
    def get_caller(metaheaders, column_header, name):
        return caller
    return get_caller

class TagTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr


    def test_build_vcf_readers(self):
        vcf_content ='''##source=strelka
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
'''
        vcf_content = vcf_content.replace('|',"\t")
        
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf",vcf_content)
            input_dir.write("B.vcf",vcf_content)
            
            vcf_readers = tag._build_vcf_readers(input_dir.path)
            
            self.assertEqual("A.vcf", vcf_readers[0].name)
            self.assertEqual(["##source=strelka"], vcf_readers[0]._metaheaders)
            self.assertEqual("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR",
                              vcf_readers[0].column_header)
            self.assertEqual("B.vcf", vcf_readers[1].name)
            self.assertEqual(["##source=strelka"], vcf_readers[1]._metaheaders)
            self.assertEqual("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR",
                              vcf_readers[1].column_header)
            self.assertEqual(2, len(vcf_readers))

    def test_build_vcf_to_caller_multipleVcfBuildsDict(self):
        vcf_content ='''##source=strelka
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
'''
        vcf_content = vcf_content.replace('|',"\t")
        
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf",vcf_content)
            input_dir.write("B.vcf",vcf_content)
            
            vcf_readers = tag._build_vcf_readers(input_dir.path)

            actual_dict = tag._build_vcf_readers_to_writers(vcf_readers, output_dir.path)
            
            actual_set = set([actual_dict.keys()[0].name,actual_dict.keys()[1].name])
            expected_set = set(["A.vcf","B.vcf"])
            self.assertTrue(actual_set.intersection(expected_set)==expected_set)
            
            actual_set = set([actual_dict.values()[0].output_filepath,actual_dict.values()[1].output_filepath])
            expected_set = set([os.path.join(output_dir.path , "B.jacquardTags.vcf"),os.path.join(output_dir.path , "B.jacquardTags.vcf")])
            self.assertTrue(actual_set.intersection(expected_set)==expected_set)


    def test_build_vcf_to_caller_multipleVcfLogs(self):
        
        vcf_content ='''##source=strelka
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
'''
        vcf_content = vcf_content.replace('|',"\t")
        
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf",vcf_content)
            input_dir.write("B.vcf",vcf_content)
            
            vcf_readers = tag._build_vcf_readers(input_dir.path)

            tag._build_vcf_readers_to_writers(vcf_readers, output_dir.path)
    
            output_lines = self.output.getvalue().rstrip().split("\n")
            self.assertEquals(3, len(output_lines))
            line_iter = iter(output_lines)
            self.assertEquals("DEBUG: VCF [A.vcf] recognized by caller [Strelka]", line_iter.next())
            self.assertEquals("DEBUG: VCF [B.vcf] recognized by caller [Strelka]", line_iter.next())
            self.assertEquals("INFO: Recognized [2] Strelka file(s)", line_iter.next())



    def test_build_vcf_readers_exceptionIsRaisedDetailsLogged(self):
        vcf_content ='''##foo
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
'''
        vcf_content = vcf_content.replace('|',"\t")
        
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf",vcf_content)
            input_dir.write("B.vcf",vcf_content)

            self.assertRaisesRegexp(jacquard_utils.JQException,
                                    "VCF files could not be parsed",
                                    tag._build_vcf_readers,
                                    input_dir.path)

    def test_tag_files(self):
        reader = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader")
        reader.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        writer = MockWriter()
        vcf_readers_to_writers = {reader: writer}
        execution_context = []
        tag.tag_files(vcf_readers_to_writers, execution_context)
        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertEquals(["##originalMeta1", 
                           "##originalMeta2",
                           "##jacquard.tag.caller=MockCaller",
                           "##mockCallerMetaheader1",
                           "#columnHeader",
                           "foo"], writer.lines())
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)
        
    def test_execute(self):
        vcf_content1 = '''##source=strelka
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
'''
        vcf_content2 = '''##source=VarScan2
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
'''
        vcf_content1 = vcf_content1.replace('|',"\t") 
        vcf_content2 = vcf_content2.replace('|',"\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf",vcf_content1)
            input_dir.write("B.vcf",vcf_content2)

            args = Namespace(input_dir=input_dir.path, 
                             output_dir=output_dir.path)

            tag.execute(args,[])

            output_dir.check("A.jacquardTags.vcf", "B.jacquardTags.vcf")
            
            file_content1 = output_dir.read("A.jacquardTags.vcf")
            file_content2 = output_dir.read("B.jacquardTags.vcf")
            
            print file_content1
            print file_content2
            
            #TODO cgates: fix somatic tags format to work
            self.assertTrue(0==1)
            
        
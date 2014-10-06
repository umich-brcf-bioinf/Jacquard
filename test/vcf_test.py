# pylint: disable=C0103,C0301,R0903,R0904
import os
import unittest
import sys
from StringIO import StringIO
from testfixtures import TempDirectory
from jacquard.vcf import VcfRecord, VcfReader, VcfWriter


class MockCaller(object):
    def __init__(self, name="MockCaller"):
        self.name = name


def build_mock_get_caller_method(caller):
    def get_caller(metaheaders, column_header):
        return caller
    return get_caller


class VcfRecordTestCase(unittest.TestCase):
    def testInit(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SAMPLE1|SAMPLE2\n".replace('|',"\t")
        record = VcfRecord(input_line)
        self.assertEquals("CHROM", record.chrom)
        self.assertEquals("POS", record.pos)
        self.assertEquals("ID", record.id)
        self.assertEquals("REF", record.ref)
        self.assertEquals("ALT", record.alt)
        self.assertEquals("QUAL", record.qual)
        self.assertEquals("FILTER", record.filter)
        self.assertEquals("INFO", record.info)
        self.assertEquals("FORMAT", record.format)
        self.assertEquals(["SAMPLE1","SAMPLE2"], record.samples)
        
    def testInit_format_set(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord(input_line)
        self.assertEquals(["F1","F2","F3"], record.format_set)

    def testInit_sample_dict(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord(input_line)
        self.assertEquals([0, 1], record.sample_dict.keys())
        self.assertEquals({"F1":"SA.1","F2":"SA.2","F3":"SA.3"}, record.sample_dict[0])
        self.assertEquals({"F1":"SB.1","F2":"SB.2","F3":"SB.3"}, record.sample_dict[1])

    def testInsertField(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord(input_line)
        record.insert_format_field("inserted",{0:0.6,1:0.5})
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:inserted|SA.1:SA.2:SA.3:0.6|SB.1:SB.2:SB.3:0.5\n".replace('|',"\t")
        self.assertEquals(expected, record.asText())

    def testInsertField_failsOnInvalidSampleDict(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord(input_line)
        self.assertRaises(KeyError, record.insert_format_field, "inserted", {0:0.6})
        self.assertRaises(KeyError, record.insert_format_field, "inserted", {0:0.6, 3:0.6})
        self.assertRaises(KeyError, record.insert_format_field, "inserted", {0:0.6, 1:0.6, 42:0.6})

    def testInsertField_failsOnExistingField(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord(input_line)
        self.assertRaises(KeyError, record.insert_format_field, "F1", {0:0.6, 1:0.6})


class VcfReaderTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr        
    
    def test_init(self):
        mock_caller = MockCaller()
        get_caller = build_mock_get_caller_method(mock_caller)

        with TempDirectory() as input_dir:
            input_dir.write("A.vcf","##source=strelka\n##foobarbaz\n#CHROM\tNORMAL\tTUMOR\n123\n456\n")
            file_path = os.path.join(input_dir.path, "A.vcf")
            reader = VcfReader(file_path, get_caller)

            self.assertEquals(file_path, reader.file_path)
            self.assertEquals("A.vcf", reader.name)
            self.assertEquals("#CHROM\tNORMAL\tTUMOR", reader.column_header)
            self.assertEquals(["##source=strelka", "##foobarbaz"], reader.metaheaders)
            self.assertEquals(mock_caller, reader.caller)

    def test_metaheadersAreImmutable(self):
        get_caller = build_mock_get_caller_method(MockCaller())
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf","##source=strelka\n##foobarbaz\n#CHROM\tNORMAL\tTUMOR\n123\n456\n")
            file_name = os.path.join(input_dir.path, "A.vcf")
            reader = VcfReader(file_name, get_caller)
            
            original_length = len(reader.metaheaders)
            reader.metaheaders.append("foo")
            
            self.assertEquals(original_length, len(reader.metaheaders))

    def test_initLogsVcfCallers(self):
        get_caller = build_mock_get_caller_method(MockCaller())
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf","##source=strelka\n##foobarbaz\n#CHROM\tNORMAL\tTUMOR\n123\n456\n")
            file_name = os.path.join(input_dir.path, "A.vcf")
            VcfReader(file_name, get_caller)

        actual_output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertEquals(1, len(actual_output_lines))
        self.assertRegexpMatches(actual_output_lines[0], "DEBUG.*A\.vcf.*MockCaller")

    def test_vcf_records(self):
        get_caller = build_mock_get_caller_method(MockCaller())
        vcf_content ='''##source=strelka
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
'''
        vcf_content = vcf_content.replace('|',"\t")
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf", vcf_content)
            file_name = os.path.join(input_dir.path, "A.vcf")
            reader = VcfReader(file_name, get_caller)

            actual_vcf_records = []
            reader.open()
            for vcf_record in reader.vcf_records():
                actual_vcf_records.append(vcf_record)
            reader.close()

        self.assertEquals(2, len(actual_vcf_records))
        self.assertEquals('chr1', actual_vcf_records[0].chrom)
        self.assertEquals('chr2', actual_vcf_records[1].chrom)


class VcfWriterTestCase(unittest.TestCase):
    def test_write(self):
        with TempDirectory() as output_dir:
            file_path = os.path.join(output_dir.path,"test.tmp")

            writer = VcfWriter(file_path)
            writer.open()
            writer.write("A")
            writer.write("B\n")
            writer.write("CD\n")
            writer.close()

            actual_output = output_dir.read('test.tmp')
            expected_output = "AB|CD|".replace('|', os.linesep)
            self.assertEquals(expected_output, actual_output)

            


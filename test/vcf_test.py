# pylint: disable=C0103,C0301,R0903,R0904
import os
import unittest
import sys
from StringIO import StringIO
from testfixtures import TempDirectory
from jacquard.vcf import VcfRecord, VcfReader, FileWriter, FileReader
import jacquard.utils as utils


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

class VcfRecordTestCase(unittest.TestCase):
    def testInit(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SAMPLE1|SAMPLE2\n".replace('|',"\t")
        record = VcfRecord.parse_record(input_line)
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
        record = VcfRecord.parse_record(input_line)
        self.assertEquals(["F1","F2","F3"], record.format_set)

    def testInit_sample_dict(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord.parse_record(input_line)
        self.assertEquals([0, 1], record.sample_dict.keys())
        self.assertEquals({"F1":"SA.1","F2":"SA.2","F3":"SA.3"}, record.sample_dict[0])
        self.assertEquals({"F1":"SB.1","F2":"SB.2","F3":"SB.3"}, record.sample_dict[1])

    def testInsertField(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord.parse_record(input_line)
        record.insert_format_field("inserted",{0:0.6,1:0.5})
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:inserted|SA.1:SA.2:SA.3:0.6|SB.1:SB.2:SB.3:0.5\n".replace('|',"\t")
        self.assertEquals(expected, record.asText())

    def testInsertField_failsOnInvalidSampleDict(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord.parse_record(input_line)
        self.assertRaises(KeyError, record.insert_format_field, "inserted", {0:0.6})
        self.assertRaises(KeyError, record.insert_format_field, "inserted", {0:0.6, 3:0.6})
        self.assertRaises(KeyError, record.insert_format_field, "inserted", {0:0.6, 1:0.6, 42:0.6})

    def testInsertField_failsOnExistingField(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord.parse_record(input_line)
        self.assertRaises(KeyError, record.insert_format_field, "F1", {0:0.6, 1:0.6})

    def testCreateKey(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        self.assertEquals("CHROM_POS_REF_ALT",VcfRecord.parse_record(input_line).key)
    
    def testInfoDict(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|k1=v1;k2=v2|F|S\n".replace('|',"\t")
        vcf_record = VcfRecord.parse_record(input_line)
        self.assertEquals({"k1":"v1","k2":"v2"},vcf_record.get_info_dict())
        
    def testEQ(self):
        base = VcfRecord.parse_record("A|1|ID|C|D|QUAL|FILTER|INFO|F|S\n".replace('|',"\t"))
        base_equivalent = VcfRecord.parse_record("A|1|ID|C|D|QUAL|FILTER||foo|S\n".replace('|',"\t"))
        self.assertEquals(base, base_equivalent)
        different_chrom = VcfRecord.parse_record("Z|1|ID|C|D|QUAL|FILTER||foo|S\n".replace('|',"\t"))
        self.assertNotEquals(base, different_chrom)
        different_pos = VcfRecord.parse_record("A|2|ID|C|D|QUAL|FILTER||foo|S\n".replace('|',"\t"))
        self.assertNotEquals(base, different_pos)
        different_ref = VcfRecord.parse_record("A|1|ID|Z|D|QUAL|FILTER||foo|S\n".replace('|',"\t"))
        self.assertNotEquals(base, different_ref)        
        different_alt = VcfRecord.parse_record("A|1|ID|C|Z|QUAL|FILTER||foo|S\n".replace('|',"\t"))
        self.assertNotEquals(base, different_alt)
        
    def testHash(self):
        base = VcfRecord.parse_record("A|B|ID|C|D|QUAL|FILTER|INFO|F|S\n".replace('|',"\t"))
        base_equivalent = VcfRecord.parse_record("A|B|ID|C|D|QUAL|FILTER||foo|S\n".replace('|',"\t"))
        self.assertEquals(base.__hash__(), base_equivalent.__hash__())
        record_set = set()
        record_set.add(base)
        record_set.add(base_equivalent)
        self.assertEquals(1,len(record_set))
        
    def testCompare(self):
        expected_records = [VcfRecord.parse_record("1|1|ID|A|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("1|1|ID|A|A|QUAL|FILTER||foo|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("1|1|ID|A|C|QUAL|FILTER|INFO|F|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("1|1|ID|C|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("1|2|ID|A|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("2|1|ID|A|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t"))]
        
        input_records = expected_records[::-1]
        
        self.assertEquals(expected_records, sorted(input_records))
        
    def testCompare_orderingByNumericChromAndPos(self):
        expected_records = [VcfRecord.parse_record("1|1|ID|A|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("2|1|ID|A|A|QUAL|FILTER||foo|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("10|1|ID|A|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("11|1|ID|C|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("20|1|ID|A|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("M|1|ID|A|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("X|1|ID|A|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t"))]
        
        input_records = expected_records[::-1]
        
        self.assertEquals(expected_records, sorted(input_records))
        
    def testCompare_nonNumericChrom(self):
        expected_records = [VcfRecord.parse_record("chr2|1|ID|A|A|QUAL|FILTER|INFO|F|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("chr5|1|ID|A|A|QUAL|FILTER||foo|S\n".replace('|',"\t")),
                   VcfRecord.parse_record("10|1|ID|A|C|QUAL|FILTER|INFO|F|S\n".replace('|',"\t"))]
        
        input_records = expected_records[::-1]
        
        self.assertEquals(expected_records, sorted(input_records))
     
    def test_empty_record(self):
        base = VcfRecord.parse_record("chr2|1|ID|A|C|QUAL|FILTER|INFO|F|S\n".replace('|',"\t"))

        empty_record = base.get_empty_record()
        
        expected_record = VcfRecord(chrom="chr2",pos="1", ref="A", alt="C")
        self.assertEquals(expected_record.asText(), empty_record.asText())
               
class VcfReaderTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr        
    
    
    def test_init(self):
        file_contents = ["##metaheader1\n",
                         "##metaheader2\n",
                         "#columnHeader\n",
                         "record1\n",
                         "record2"]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)

        actual_vcf_reader = VcfReader(mock_reader)
        
        self.assertEquals("my_dir/my_file.txt", actual_vcf_reader.input_filepath)
        self.assertEquals("my_file.txt", actual_vcf_reader.file_name)
        self.assertEquals("#columnHeader", actual_vcf_reader.column_header)
        self.assertEquals(["##metaheader1", "##metaheader2"], actual_vcf_reader.metaheaders)

    def test_metaheadersAreImmutable(self):
        file_contents = ["##metaheader1\n",
                         "##metaheader2\n",
                         "#columnHeader\n",
                         "record1\n",
                         "record2"]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)
        
        original_length = len(reader.metaheaders)
        reader.metaheaders.append("foo")
        
        self.assertEquals(original_length, len(reader.metaheaders))

    def test_vcf_records(self):
        file_contents = ["##metaheader1\n",
                         "##metaheader2\n",
                         "#columnHeader\n",
                         "chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR\n".replace("|","\t"),
                         "chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR".replace("|","\t")]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)

        actual_vcf_records = []
        reader.open()
        for vcf_record in reader.vcf_records():
            actual_vcf_records.append(vcf_record)
        reader.close()

        self.assertEquals(2, len(actual_vcf_records))
        self.assertEquals('chr1', actual_vcf_records[0].chrom)
        self.assertEquals('chr2', actual_vcf_records[1].chrom)
        self.assertTrue(mock_reader.open_was_called)
        self.assertTrue(mock_reader.close_was_called)

        
    def test_noColumnHeaders(self):
        mock_reader = MockFileReader("my_dir/my_file.txt", ["##metaheader\n"])
        self.assertRaises(utils.JQException, VcfReader, mock_reader)
            
    def test_noMetaheaders(self):
        mock_reader = MockFileReader("my_dir/my_file.txt", ["#columnHeader\n"])
        self.assertRaises(utils.JQException, VcfReader, mock_reader)
            
        
class VcfWriterTestCase(unittest.TestCase):
    def test_write(self):
        with TempDirectory() as output_dir:
            file_path = os.path.join(output_dir.path,"test.tmp")

            writer = FileWriter(file_path)
            writer.open()
            writer.write("A")
            writer.write("B\n")
            writer.write("CD\n")
            writer.close()

            actual_output = output_dir.read('test.tmp')
            expected_output = "AB|CD|".replace('|', os.linesep)
            self.assertEquals(expected_output, actual_output)
            
class FileReaderTestCase(unittest.TestCase):
    def test_equality(self):
        self.assertEquals(FileReader("foo"),FileReader("foo"))
        self.assertNotEquals(FileReader("foo"),FileReader("bar"))
        self.assertNotEquals(FileReader("foo"), 1)

    def test_hashable(self):
        s = set([FileReader("foo")])
        s.add(FileReader("foo"))
        self.assertEquals(1, len(s))
    
    def test_read_lines(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.tmp","1\n2\n3")
            reader = FileReader(os.path.join(input_dir.path, "A.tmp"))
            reader.open()
            actual_lines = [line for line in reader.read_lines()]
            reader.close()

            self.assertEquals(["1\n","2\n","3"], actual_lines)

class FileWriterTestCase(unittest.TestCase):
    def test_equality(self):
        self.assertEquals(FileWriter("foo"),FileWriter("foo"))
        self.assertNotEquals(FileWriter("foo"),FileWriter("bar"))
        self.assertNotEquals(FileWriter("foo"), 1)
            
    def test_hashable(self):
        s = set([FileWriter("foo")])
        s.add(FileWriter("foo"))
        self.assertEquals(1, len(s))
            
    def test_write_lines(self):
        with TempDirectory() as output_dir:
            writer = FileWriter(os.path.join(output_dir.path, "A.tmp"))
            writer.open()
            writer.write("1\n2\n")
            writer.write("3")      
            writer.close()

            actual_file = open(os.path.join(output_dir.path, "A.tmp"))
            actual_output = actual_file.readlines()
            actual_file.close()
            
            self.assertEquals(["1\n","2\n","3"], actual_output)
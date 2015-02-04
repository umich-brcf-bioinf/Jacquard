#pylint: disable=line-too-long,too-many-public-methods,invalid-name
#pylint: disable=missing-docstring,protected-access,too-few-public-methods
from StringIO import StringIO
from collections import OrderedDict
from jacquard.vcf import VcfRecord, VcfReader, FileWriter, FileReader
from testfixtures import TempDirectory
import jacquard.utils as utils
import os
import sys
import test.test_case as test_case
import unittest

class MockFileReader(object):
    def __init__(self, input_filepath="/foo/mockFileReader.txt", content=None):
        self.input_filepath = input_filepath
        self.file_name = os.path.basename(input_filepath)
        if content is None:
            self.content = []
        else:
            self.content = content
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

#TODO: (cgates) Fix tests to that so that they do not rely on parse_record and asText.
class VcfRecordTestCase(test_case.JacquardBaseTestCase):
    def test_parse_record(self):
        sample_names = ["SampleA", "SampleB"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FOO:BAR|SA_foo:SA_bar|SB_foo:SB_bar\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        self.assertEquals("CHROM", record.chrom)
        self.assertEquals("POS", record.pos)
        self.assertEquals("ID", record.id)
        self.assertEquals("REF", record.ref)
        self.assertEquals("ALT", record.alt)
        self.assertEquals("QUAL", record.qual)
        self.assertEquals("FILTER", record.filter)
        self.assertEquals("INFO", record.info)

    def test_format_tags(self):
        sample_names = ["SampleA", "SampleB"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        self.assertEquals(set(["F1", "F2", "F3"]), record.format_tags)

    def test_format_tags_emptyWhenNoSamples(self):
        sample_names = []
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        self.assertEquals(set(), record.format_tags)

    def test_format_field(self):
        sample_names = ["SA", "SB"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F3:F1:F2|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        self.assertEquals("F3:F1:F2", record._format_field())

    def test_format_field_emptyWhenNoSamples(self):
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO\n")
        record = VcfRecord.parse_record(input_line, [])
        self.assertEquals(".", record._format_field())

    def test_format_field_preservesOrderWhenAddingNewTags(self):
        sample_names = ["SA", "SB"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F3:F1:F2|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        record.add_sample_tag_value("Z4", {"SA" : "SA.4", "SB" : "SB.4"})
        record.add_sample_tag_value("A5", {"SA"  :"SA.A5", "SB" : "SB.A5"})
        self.assertEquals("F3:F1:F2:Z4:A5", record._format_field())

    def test_parse_record_sample_dict(self):
        sample_names = ["SampleA", "SampleB"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        self.assertEquals(["SampleA", "SampleB"], record.sample_tag_values.keys())
        self.assertEquals({"F1":"SA.1", "F2":"SA.2", "F3":"SA.3"}, record.sample_tag_values["SampleA"])
        self.assertEquals({"F1":"SB.1", "F2":"SB.2", "F3":"SB.3"}, record.sample_tag_values["SampleB"])

    def test_sample_tag_values(self):
        sample_tag_values = VcfRecord._sample_tag_values(["sampleA", "sampleB"],
                                                         "foo:bar",
                                                         ["SA_foo:SA_bar", "SB_foo:SB_bar"])
        self.assertEquals({"foo":"SA_foo", "bar":"SA_bar"}, sample_tag_values["sampleA"])
        self.assertEquals({"foo":"SB_foo", "bar":"SB_bar"}, sample_tag_values["sampleB"])

    def test_filter_sample_tag_values(self):
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FOO:BAR|A:B|C:D\n")
        record = VcfRecord.parse_record(input_line, sample_names=["sampleA", "sampleB"])
        desired_tags = ["FO*"]
        record.filter_sample_tag_values(desired_tags)
        self.assertEquals({"sampleA":{"FOO":"A"}, "sampleB":{"FOO":"C"}}, record.sample_tag_values)

        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FOO:BAR|A:B|C:D\n")
        record = VcfRecord.parse_record(input_line, sample_names=["sampleA", "sampleB"])
        desired_tags = ["BA*"]
        record.filter_sample_tag_values(desired_tags)
        self.assertEquals({"sampleA":{"BAR":"B"}, "sampleB":{"BAR":"D"}}, record.sample_tag_values)

    def test_sample_tag_values_emptyDictWhenExplicitNullSampleData(self):
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|.|.|.\n")
        record = VcfRecord.parse_record(input_line, sample_names=["sampleA", "sampleB"])
        self.assertEquals(["sampleA", "sampleB"], record.sample_tag_values.keys())
        self.assertEquals({}, record.sample_tag_values["sampleA"])
        self.assertEquals({}, record.sample_tag_values["sampleB"])

    def test_sample_tag_values_whenSparseSampleData(self):
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FOO|.|.\n")
        record = VcfRecord.parse_record(input_line, sample_names=["sampleA", "sampleB"])
        self.assertEquals(["sampleA", "sampleB"], record.sample_tag_values.keys())
        self.assertEquals(OrderedDict({"FOO":"."}), record.sample_tag_values["sampleA"])
        self.assertEquals(OrderedDict({"FOO":"."}), record.sample_tag_values["sampleB"])

    def test_sample_tag_values_emptyDictWhenNoSampleData(self):
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|||\n")
        record = VcfRecord.parse_record(input_line, sample_names=["sampleA", "sampleB"])
        self.assertEquals(["sampleA", "sampleB"], record.sample_tag_values.keys())
        self.assertEquals({}, record.sample_tag_values["sampleA"])
        self.assertEquals({}, record.sample_tag_values["sampleB"])

    def test_sample_tag_values_emptyDictWhenNoSamples(self):
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO\n")
        record = VcfRecord.parse_record(input_line, sample_names=["sampleA", "sampleB"])
        self.assertEquals({}, record.sample_tag_values)

    def test_parse_record_initsSampleTagValues(self):
        sample_names = ["SampleA", "SampleB"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        self.assertEquals(["SampleA", "SampleB"], record.sample_tag_values.keys())
        self.assertEquals({"F1":"SA.1", "F2":"SA.2", "F3":"SA.3"}, record.sample_tag_values["SampleA"])
        self.assertEquals({"F1":"SB.1", "F2":"SB.2", "F3":"SB.3"}, record.sample_tag_values["SampleB"])

    def test_sample_tag_values_preservesSampleOrder(self):
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|||\n")
        record = VcfRecord.parse_record(input_line, sample_names=["sampleB", "sampleA"])
        self.assertEquals(["sampleB", "sampleA"], record.sample_tag_values.keys())

    def test_add_sample_format_value(self):
        sample_names = ["SampleA", "SampleB"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        record.add_sample_tag_value("inserted", {"SampleB":"insertedValueB", "SampleA":"insertedValueA"})
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:inserted|SA.1:SA.2:SA.3:insertedValueA|SB.1:SB.2:SB.3:insertedValueB\n")
        self.assertEquals(expected, record.asText())

    def test_insert_format_field_failsOnInvalidSampleDict(self):
        sample_names = ["SampleA", "SampleB"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        self.assertRaises(KeyError, record.add_sample_tag_value, "inserted", {"SampleA":0.6})
        self.assertRaises(KeyError, record.add_sample_tag_value, "inserted", {"SampleA":0.6, "SampleZ":0.6})
        self.assertRaises(KeyError, record.add_sample_tag_value, "inserted", {"SampleA":0.6, "SampleB":0.6, "SampleZ":0.6})

    def test_insert_format_field_failsOnExistingField(self):
        sample_names = ["SampleA", "SampleB"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        self.assertRaises(KeyError, record.add_sample_tag_value, "F1", {"SampleA":0.6, "SampleB":0.6})

    def test_get_info_dict_empty(self):
        vcf_record = VcfRecord("chr1", "42", "A", "C", info = "")
        self.assertEquals({}, vcf_record.info_dict)

    def test_get_info_dict_null(self):
        vcf_record = VcfRecord("chr1", "42", "A", "C", info = ".")
        self.assertEquals({}, vcf_record.info_dict)

    def test_add_info_field_assignedField(self):
        sample_names = ["SampleA"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|k1=v1;k2=v2;baz|F|S\n")
        vcf_record = VcfRecord.parse_record(input_line, sample_names)
        vcf_record.add_info_field("foo=bar")
        self.assertEquals({"k1": "v1", "k2": "v2", "baz": "baz", "foo": "bar"}, vcf_record.info_dict)

    def test_add_info_field_nonAssignedField(self):
        sample_names = ["SampleA"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|k1=v1;k2=v2;baz|F|S\n")
        vcf_record = VcfRecord.parse_record(input_line, sample_names)
        vcf_record.add_info_field("foo")
        self.assertEquals({"k1": "v1", "k2": "v2", "baz": "baz", "foo": "foo"}, vcf_record.info_dict)

    def test_join_info_fields(self):
        sample_names = ["SampleA"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|k1=v1;k2=v2;baz|F|S\n")
        vcf_record = VcfRecord.parse_record(input_line, sample_names)
        vcf_record._join_info_fields()
        self.assertEquals("k2=v2;k1=v1;baz", vcf_record.info)

    def test_join_info_fields_nullValues(self):
        sample_names = ["SampleA"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|.|F|S\n")
        vcf_record = VcfRecord.parse_record(input_line, sample_names)
        vcf_record._join_info_fields()
        self.assertEquals(".", vcf_record.info)

        vcf_record = VcfRecord.parse_record(input_line, sample_names)
        vcf_record.add_info_field("foo")
        vcf_record._join_info_fields()
        self.assertEquals("foo", vcf_record.info)

    def test_asText(self):
        sampleA = OrderedDict({"F1":"SA.1", "F2":"SA.2", "F3":"SA.3"})
        sampleB = OrderedDict({"F1":"SB.1", "F2":"SB.2", "F3":"SB.3"})
        sample_tag_values = OrderedDict({"SampleA":sampleA, "SampleB":sampleB})
        record = VcfRecord("CHROM", "POS", "REF", "ALT", "ID", "QUAL", "FILTER", "INFO", sample_tag_values)
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        self.assertEquals(expected, record.asText())

    def test_asTextWhenEmptyFormatField(self):
        sampleA = OrderedDict({})
        sampleB = OrderedDict({})
        sample_tag_values = OrderedDict({"SampleA":sampleA, "SampleB":sampleB})
        record = VcfRecord("CHROM", "POS", "REF", "ALT", "ID", "QUAL", "FILTER", "INFO", sample_tag_values)
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|.|.|.\n")
        self.assertEquals(expected, record.asText())

    def test_equals(self):
        sample_names = ["sampleA"]
        base = VcfRecord.parse_record(self.entab("A|1|ID|C|D|QUAL|FILTER|INFO|F|S\n"), sample_names)
        base_equivalent = VcfRecord.parse_record(self.entab("A|1|ID|C|D|QUAL|FILTER||foo|S\n"), sample_names)
        self.assertEquals(base, base_equivalent)
        different_chrom = VcfRecord.parse_record(self.entab("Z|1|ID|C|D|QUAL|FILTER||foo|S\n"), sample_names)
        self.assertNotEquals(base, different_chrom)
        different_pos = VcfRecord.parse_record(self.entab("A|2|ID|C|D|QUAL|FILTER||foo|S\n"), sample_names)
        self.assertNotEquals(base, different_pos)
        different_ref = VcfRecord.parse_record(self.entab("A|1|ID|Z|D|QUAL|FILTER||foo|S\n"), sample_names)
        self.assertNotEquals(base, different_ref)
        different_alt = VcfRecord.parse_record(self.entab("A|1|ID|C|Z|QUAL|FILTER||foo|S\n"), sample_names)
        self.assertNotEquals(base, different_alt)

    def testHash(self):
        sample_names = ["sampleA"]
        base = VcfRecord.parse_record(self.entab("A|B|ID|C|D|QUAL|FILTER|INFO|F|S\n"), sample_names)
        base_equivalent = VcfRecord.parse_record(self.entab("A|B|ID|C|D|QUAL|FILTER||foo|S\n"), sample_names)
        self.assertEquals(base.__hash__(), base_equivalent.__hash__())
        record_set = set()
        record_set.add(base)
        record_set.add(base_equivalent)
        self.assertEquals(1, len(record_set))

    def testCompare(self):
        sample_names = ["SampleA"]
        expected_records = [VcfRecord.parse_record(self.entab("1|1|ID|A|A|QUAL|FILTER|INFO|F|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("1|1|ID|A|A|QUAL|FILTER||foo|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("1|1|ID|A|C|QUAL|FILTER|INFO|F|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("1|1|ID|C|A|QUAL|FILTER|INFO|F|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("1|2|ID|A|A|QUAL|FILTER|INFO|F|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("2|1|ID|A|A|QUAL|FILTER|INFO|F|S\n"), sample_names)]

        input_records = expected_records[::-1]

        self.assertEquals(expected_records, sorted(input_records))

    def testCompare_orderingByNumericChromAndPos(self):
        sample_names = ["SampleA"]
        expected_records = [VcfRecord.parse_record(self.entab("1|1|ID|A|A|QUAL|FILTER|INFO|F|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("2|1|ID|A|A|QUAL|FILTER||foo|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("10|1|ID|A|A|QUAL|FILTER|INFO|F|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("11|1|ID|C|A|QUAL|FILTER|INFO|F|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("20|1|ID|A|A|QUAL|FILTER|INFO|F|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("M|1|ID|A|A|QUAL|FILTER|INFO|F|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("X|1|ID|A|A|QUAL|FILTER|INFO|F|S\n"), sample_names)]

        input_records = expected_records[::-1]

        self.assertEquals(expected_records, sorted(input_records))

    def testCompare_nonNumericChrom(self):
        sample_names = ["SampleA"]
        expected_records = [VcfRecord.parse_record(self.entab("chr2|1|ID|A|A|QUAL|FILTER|INFO|F|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("chr5|1|ID|A|A|QUAL|FILTER||foo|S\n"), sample_names),
                            VcfRecord.parse_record(self.entab("10|1|ID|A|C|QUAL|FILTER|INFO|F|S\n"), sample_names)]

        input_records = expected_records[::-1]

        self.assertEquals(expected_records, sorted(input_records))

    def test_empty_record(self):
        sample_names = ["SampleA"]
        base = VcfRecord.parse_record(self.entab("chr2|1|ID|A|C|QUAL|FILTER|INFO|F|S\n"), sample_names)

        empty_record = base.get_empty_record()

        expected_record = VcfRecord(chrom="chr2", pos="1", ref="A", alt="C")
        self.assertEquals(expected_record.asText(), empty_record.asText())

class VcfReaderTestCase(test_case.JacquardBaseTestCase):
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
        self.assertEquals(("##metaheader1", "##metaheader2"), actual_vcf_reader.metaheaders)
        self.assertEquals([], actual_vcf_reader.sample_names)

    def test_init_sampleNamesInitialized(self):
        file_contents = ["##metaheader1\n",
                         "##metaheader2\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleA|SampleB\n"),
                         "record1\n",
                         "record2"]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)

        actual_vcf_reader = VcfReader(mock_reader)
        self.assertEquals(["SampleA", "SampleB"], actual_vcf_reader.sample_names)

    def test_format_metaheaders(self):
        file_contents = ["##metaheader1\n",
                         "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n",
                         "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)

        expected_metaheaders = {"AF" : '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
                                "DP" : '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">'}

        self.assertEquals(expected_metaheaders, reader.format_metaheaders)

    def test_info_metaheaders(self):
        file_contents = ["##metaheader1\n",
                         "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n",
                         "##INFO=<ID=SNP,Number=1,Type=Integer,Description=\"snp\">\n",
                         "##INFO=<ID=FOO,Number=1,Type=Integer,Description=\"foo\">\n",
                         "##INFO=<ID=BAR,Number=1,Type=Integer,Description=\"bar\">\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr1|1|.|A|C|.|.|SNP;BAR|FORMAT|NORMAL|TUMOR\n"),
                         self.entab("chr2|1|.|A|C|.|.|BAR|FORMAT|NORMAL|TUMOR")]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)

        expected_metaheaders = {"SNP" : "##INFO=<ID=SNP,Number=1,Type=Integer,Description=\"snp\">",
                                "FOO" : "##INFO=<ID=FOO,Number=1,Type=Integer,Description=\"foo\">",
                                "BAR" : "##INFO=<ID=BAR,Number=1,Type=Integer,Description=\"bar\">"}

        self.assertEquals(expected_metaheaders, reader.info_metaheaders)

    def test_filter_metaheaders(self):
        file_contents = ["##metaheader1\n",
                         "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n",
                         "##INFO=<ID=SNP,Number=1,Type=Integer,Description=\"snp\">\n",
                         "##FILTER=<ID=.,Number=1,Type=Integer,Description=\"foo\">\n",
                         "##FILTER=<ID=PASS,Number=1,Type=Integer,Description=\"bar\">\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)

        expected_metaheaders = {"." : "##FILTER=<ID=.,Number=1,Type=Integer,Description=\"foo\">",
                                "PASS" : "##FILTER=<ID=PASS,Number=1,Type=Integer,Description=\"bar\">",}

        self.assertEquals(expected_metaheaders, reader.filter_metaheaders)

    def test_format_tag_ids_ignoresRelatedFieldNames(self):
        file_contents = ["##metaheader1\n",
                         "##FORMAT=<UUID=DPX1,ID=DP1>\n",
                         "##FORMAT=<ID=DP2,UUID=DPX2>\n",
                         "##FORMAT=<ID=DP3,ID=DPX3>\n",
                         "##FORMAT=<NOID=DPX>\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)

        self.assertEquals(set(["DP1", "DP2", "DP3"]), set(reader.format_metaheaders.keys()))


    def test_format_tag_ids_idsAreUnique(self):
        file_contents = ["##metaheader1\n",
                         "##FORMAT=<ID=AF,Description='Allele Frequency 1'>\n",
                         "##FORMAT=<ID=AF,Description='Allele Frequency 2'>\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)

        self.assertEquals(["AF"], reader.format_metaheaders.keys())
        self.assertEquals("##FORMAT=<ID=AF,Description='Allele Frequency 2'>", reader.format_metaheaders["AF"])

    def test_format_tag_ids_emptyWhenNoFormatTags(self):
        file_contents = ["##metaheader1\n",
                         "##INFO=<ID=AF,Number=A,Type=Float,Description='Allele Frequency'>\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)

        self.assertEquals(0, len(reader.format_metaheaders))

    def test_format_tag_ids_immutable(self):
        file_contents = ["##metaheader1\n",
                         "##FORMAT=<ID=DP,Number=1,Type=Integer,Description='Read Depth'>\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)

        self.assertEquals(["DP"], reader.format_metaheaders.keys())
        del reader.format_metaheaders["DP"]
        self.assertEquals(["DP"], reader.format_metaheaders.keys())

    def test_vcf_records(self):
        file_contents = ["##metaheader1\n",
                         "##metaheader2\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
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

    def test_get_format_tag_list(self):
        file_contents = ['##FORMAT=<ID=GT,Number=1>\n',
                         '##FORMAT=<ID=GQ,Number=1,Description="bar">\n',
                         '#columnHeader\n',
                         'record1\n',
                         'record2']
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)

        vcf_reader = VcfReader(mock_file_reader)
        actual_format_set = vcf_reader.format_metaheaders
        expected_format_set = ["GT", "GQ"]

        self.assertEquals(expected_format_set, actual_format_set.keys())

    def test_get_info_field_list(self):
        file_contents = ['##INFO=<ID=AF,Number=1>\n',
                         '##FORMAT=<ID=GQ,Number=1,Description="bar">\n',
                         '##INFO=<ID=AA,Number=1>\n',
                         '#columnHeader\n',
                         'record1\n',
                         'record2']
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)

        vcf_reader = VcfReader(mock_file_reader)
        actual_format_set = vcf_reader.info_metaheaders
        expected_format_set = ["AA", "AF"]

        self.assertEquals(expected_format_set, actual_format_set.keys())

class VcfWriterTestCase(unittest.TestCase):
    def test_write(self):
        with TempDirectory() as output_file:
            file_path = os.path.join(output_file.path, "test.tmp")

            writer = FileWriter(file_path)
            writer.open()
            writer.write("A")
            writer.write("B\n")
            writer.write("CD\n")
            writer.close()

            actual_output = output_file.read('test.tmp')
            expected_output = "AB|CD|".replace('|', os.linesep)
            self.assertEquals(expected_output, actual_output)

class FileReaderTestCase(unittest.TestCase):
    def test_equality(self):
        self.assertEquals(FileReader("foo"), FileReader("foo"))
        self.assertNotEquals(FileReader("foo"), FileReader("bar"))
        self.assertNotEquals(FileReader("foo"), 1)

    def test_hashable(self):
        s = set([FileReader("foo")])
        s.add(FileReader("foo"))
        self.assertEquals(1, len(s))

    def test_read_lines(self):
        with TempDirectory() as input_file:
            input_file.write("A.tmp", "1\n2\n3")
            reader = FileReader(os.path.join(input_file.path, "A.tmp"))
            reader.open()
            actual_lines = [line for line in reader.read_lines()]
            reader.close()

            self.assertEquals(["1\n", "2\n", "3"], actual_lines)

class FileWriterTestCase(unittest.TestCase):
    def test_equality(self):
        self.assertEquals(FileWriter("foo"), FileWriter("foo"))
        self.assertNotEquals(FileWriter("foo"), FileWriter("bar"))
        self.assertNotEquals(FileWriter("foo"), 1)

    def test_file_name(self):
        writer = FileWriter("foo/bar/baz.tmp")
        self.assertEquals("baz.tmp", writer.file_name)

    def test_hashable(self):
        s = set([FileWriter("foo")])
        s.add(FileWriter("foo"))
        self.assertEquals(1, len(s))

    def test_write_lines(self):
        with TempDirectory() as output_file:
            writer = FileWriter(os.path.join(output_file.path, "A.tmp"))
            writer.open()
            writer.write("1\n2\n")
            writer.write("3")
            writer.close()

            actual_file = open(os.path.join(output_file.path, "A.tmp"))
            actual_output = actual_file.readlines()
            actual_file.close()

            self.assertEquals(["1\n", "2\n", "3"], actual_output)

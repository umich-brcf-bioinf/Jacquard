#pylint: disable=line-too-long,too-many-public-methods,invalid-name
#pylint: disable=missing-docstring,protected-access,too-few-public-methods
#pylint: disable=too-many-arguments,too-many-instance-attributes
from StringIO import StringIO
from collections import OrderedDict
import os
import re
import sys
import unittest

from testfixtures import TempDirectory

import jacquard.utils as utils
from jacquard.vcf import VcfRecord, VcfReader, FileWriter, FileReader
import test.test_case as test_case

class MockFileWriter(object):
    def __init__(self):
        self._content = []
        self.opened = False
        self.closed = False

    def open(self):
        self.opened = True

    def write(self, content):
        if content == None:
            return
        self._content.extend(content.splitlines())

    def lines(self):
        return self._content

    def close(self):
        self.closed = True


class MockFileReader(object):
    def __init__(self, input_filepath="/foo/mockFileReader.txt", content=None):
        self.input_filepath = input_filepath
        self.file_name = os.path.basename(input_filepath)
        if content is None:
            self.content = []
        else:
            self.content = content
        self.open_was_called = False
        self.close_was_called = False
        self.lines_to_iterate = None

    def open(self):
        self.open_was_called = True
        self.lines_to_iterate = list(self.content)

    def read_lines(self):
        for line in self.lines_to_iterate:
            yield line

    def close(self):
        self.close_was_called = True
        self.lines_to_iterate = None

class MockWriter(object):
    def __init__(self):
        self._content = []
        self.output_filepath = "foo"
        self.opened = False
        self.closed = False

    def open(self):
        self.opened = True

    def write(self, content):
        self._content.extend(content.splitlines())

    def lines(self):
        return self._content

    def close(self):
        self.closed = True

class MockVcfReader(object):
    def __init__(self,
                 input_filepath="vcfName",
                 metaheaders=None,
                 column_header='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR',
                 content=None,
                 records=None,
                 sample_names=None):

        if content is None:
            self.content = ["foo"]
        else:
            self.content = content

        if metaheaders is None:
            self.metaheaders = ["##metaheaders"]
        else:
            self.metaheaders = metaheaders

        if records:
            self.records = records
        elif content:
            self.records = [MockVcfRecord.parse_record(line) for line in self.content]
        else:
            self.records = []

        if sample_names is None:
            self.sample_names = []
        else:
            self.sample_names = sample_names
        self.file_name = input_filepath
        self.input_filepath = input_filepath
        self.column_header = column_header
        self.split_column_header = self.column_header.strip("#").split("\t")
        self.opened = False
        self.closed = False


    def open(self):
        self.opened = True

    def vcf_records(self):
        for record in self.records:
            yield record

    def _get_tag_metaheaders(self, regex_exp):
        tag_dict = {}
        for metaheader in self.metaheaders:
            tag = re.match(regex_exp, metaheader)
            if tag:
                tag_key = tag.group(1)
                tag_dict[tag_key] = metaheader.strip()

        return tag_dict

    @property
    def format_metaheaders(self):
        return dict(self._get_tag_metaheaders("^##FORMAT=.*?[<,]ID=([^,>]*)"))

    @property
    def info_metaheaders(self):
        return dict(self._get_tag_metaheaders("^##INFO=.*?[<,]ID=([^,>]*)"))

    @property
    def filter_metaheaders(self):
        return dict(self._get_tag_metaheaders("^##FILTER=.*?[<,]ID=([^,>]*)"))

    @property
    def contig_metaheaders(self):
        return dict(self._get_tag_metaheaders("^##contig=.*?[<,]ID=([^,>]*)"))

    @property
    def non_format_metaheaders(self):
        return self.metaheaders

    def close(self):
        self.closed = True


class MockCaller(object):
    def __init__(self, name="MockCaller", metaheaders=None, claimable=None):
        if claimable:
            self.claimable = claimable
        else:
            self.claimable = set()
        self.name = name
        if metaheaders:
            self.metaheaders = metaheaders
        else:
            self.metaheaders = ["##mockMetaheader1"]
        self.file_name_search = "snps|indels"

    @staticmethod
    def add_tags(vcfRecord):
        return vcfRecord

    @staticmethod
    def decorate_files(filenames, dummy):
        return filenames[0]+"foo"

    def get_new_metaheaders(self):
        return self.metaheaders

    def claim(self, file_readers):
        unclaimed = list(set(file_readers).difference(self.claimable))
        claimed = list(set(file_readers).intersection(self.claimable))
        return unclaimed, claimed

class MockVcfRecord(object):
    @classmethod
    def parse_record(cls, vcf_line):
        vcf_fields = vcf_line.rstrip().split("\t")
        chrom, pos, rid, ref, alt, qual, rfilter, info, rformat \
                = vcf_fields[0:9]
        samples = vcf_fields[9:]
        return MockVcfRecord(chrom, pos, ref, alt, rid, qual, rfilter, info,
                             rformat, samples)

    def __init__(self, chrom, pos, ref, alt,
                 vcf_id=".", qual=".", vcf_filter=".", info=".", vcf_format=".",
                 samples=None):
        self.chrom = chrom
        self.pos = pos
        self.id = vcf_id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = vcf_filter
        self.info = info
        self.format = vcf_format
        if samples is None:
            self.samples = []
        else:
            self.samples = samples

        tags = self.format.split(":")
        self.format_set = tags

        self.sample_dict = {}
        for i, sample in enumerate(self.samples):
            values = sample.split(":")
            self.sample_dict[i] = OrderedDict(zip(tags, values))

    def get_empty_record(self):
        return MockVcfRecord(self.chrom, self.pos, self.ref, self.alt)

    def get_info_dict(self):
        info_dict = {}

        for key_value in self.info.split(";"):
            if "=" in key_value:
                key, value = key_value.split("=")
                info_dict[key] = value
            else:
                info_dict[key_value] = key_value

        return info_dict

    def text(self):
        stringifier = [self.chrom, self.pos, self.id, self.ref, self.alt,
                       self.qual, self.filter, self.info,
                       ":".join(self.format_set)]

        for key in self.sample_dict:
            stringifier.append(":".join(self.sample_dict[key].values()))

        return "\t".join(stringifier) + "\n"

    def __eq__(self, other):
        return ("^".join([self.chrom,
                          self.pos,
                          self.ref,
                          self.alt]) == other)

class MockTag(object):
    def __init__(self, field_name=None, sample_values=None, metaheader=None):
        self.field_name = field_name
        if sample_values:
            self.sample_values = sample_values
        else:
            self.sample_values = {}
        if metaheader:
            self.metaheader = metaheader
        else:
            self.metaheader = []

    def add_tag_values(self, vcf_record):
        vcf_record.add_sample_tag_value(self.field_name, self.sample_values)

#TODO: (cgates) Fix tests to not use parse_record() and text().
class VcfRecordTestCase(test_case.JacquardBaseTestCase):
    def test_parse_record(self):
        sample_names = ["SampleA", "SampleB"]
        input_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FOO:BAR|SA_foo:SA_bar|SB_foo:SB_bar\n")
        record = VcfRecord.parse_record(input_line, sample_names)
        self.assertEquals("CHROM", record.chrom)
        self.assertEquals("POS", record.pos)
        self.assertEquals("ID", record.vcf_id)
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
        self.assertEquals(expected, record.text())

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
        vcf_record = VcfRecord("chr1", "42", "A", "C", info="")
        self.assertEquals({}, vcf_record.info_dict)

    def test_get_info_dict_null(self):
        vcf_record = VcfRecord("chr1", "42", "A", "C", info=".")
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

    def test_text(self):
        sampleA = OrderedDict({"F1":"SA.1", "F2":"SA.2", "F3":"SA.3"})
        sampleB = OrderedDict({"F1":"SB.1", "F2":"SB.2", "F3":"SB.3"})
        sample_tag_values = OrderedDict({"SampleA":sampleA, "SampleB":sampleB})
        record = VcfRecord("CHROM", "POS", "REF", "ALT", "ID", "QUAL", "FILTER", "INFO", sample_tag_values)
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n")
        self.assertEquals(expected, record.text())

    def test_asTextWhenEmptyFormatField(self):
        sampleA = OrderedDict({})
        sampleB = OrderedDict({})
        sample_tag_values = OrderedDict({"SampleA":sampleA, "SampleB":sampleB})
        record = VcfRecord("CHROM", "POS", "REF", "ALT", "ID", "QUAL", "FILTER", "INFO", sample_tag_values)
        expected = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|.|.|.\n")
        self.assertEquals(expected, record.text())

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
        self.assertEquals(expected_record.text(), empty_record.text())

    def test_add_or_replace_filter_filterReplacesPassFilter(self):
        record = VcfRecord("chr1", "42", "X", "C", vcf_filter="PASS")
        record.add_or_replace_filter("JQ_EXCLUDE")
        self.assertEquals("JQ_EXCLUDE", record.filter)

    def test_add_or_replace_filter_filterReplacesNullFilter(self):
        record = VcfRecord("chr1", "42", "X", "C", vcf_filter=".")
        record.add_or_replace_filter("JQ_EXCLUDE")
        self.assertEquals("JQ_EXCLUDE", record.filter)

    def test_add_or_replace_filter_filterReplacesEmptyFilter(self):
        record = VcfRecord("chr1", "42", "X", "C", vcf_filter="")
        record.add_or_replace_filter("JQ_EXCLUDE")
        self.assertEquals("JQ_EXCLUDE", record.filter)

    def test_add_or_replace_filter_filterAppendsFailedFilter(self):
        record = VcfRecord("chr1", "42", "XYZ", "C", vcf_filter="indelError")
        record.add_or_replace_filter("JQ_EXCLUDE")
        self.assertEquals("indelError;JQ_EXCLUDE", record.filter)

    def test_add_or_replace_filter_duplicateFilterNotAdded(self):
        record = VcfRecord("chr1", "42", "XYZ", "C", vcf_filter="JQ_EXCLUDE")
        record.add_or_replace_filter("JQ_EXCLUDE")
        self.assertEquals("JQ_EXCLUDE", record.filter)

    def test_add_or_replace_filter_filtersOnlyAppendsUnique(self):
        record = VcfRecord("chr1", "42", "XYZ", "C", vcf_filter="indelError")
        record.add_or_replace_filter("JQ_EXCLUDE")
        record.add_or_replace_filter("JQ_EXCLUDE")
        self.assertEquals("indelError;JQ_EXCLUDE", record.filter)


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

    def test_vcf_records_raisesStopIterationWhenExhausted(self):
        file_contents = ["##metaheader1\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR\n")]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)

        reader.open()
        record_iter = reader.vcf_records()
        record_iter.next()
        self.assertRaises(StopIteration,
                          record_iter.next)

    def test_vcf_records_raisesTypeErrorWhenClosed(self):
        file_contents = ["##metaheader1\n",
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR\n")]
        mock_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        reader = VcfReader(mock_reader)

        record_iter = reader.vcf_records()
        self.assertRaises(TypeError,
                          record_iter.next)

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
    def testCompare(self):
        expected_readers = [FileReader("1A.txt"),
                            FileReader("1B.txt"),
                            FileReader("2A.txt"),
                            FileReader("10A.txt"),
                            FileReader("10B.txt"),
                            FileReader("11A.txt"),
                            FileReader("11B.txt"),
                            FileReader("20A.txt"),
                            FileReader("100A.txt")]
        input_readers = expected_readers[::-1]

        self.assertEquals(expected_readers, sorted(input_readers))

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

    def test_read_lines_raisesTypeErrorWhenClosed(self):
        with TempDirectory() as input_file:
            input_file.write("A.tmp", "1\n2\n3")
            reader = FileReader(os.path.join(input_file.path, "A.tmp"))
            line_iter = reader.read_lines()
            self.assertRaises(TypeError, line_iter.next)

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

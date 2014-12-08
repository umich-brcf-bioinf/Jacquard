# pylint: disable=C0103,C0301,R0903,R0904,C0111
from argparse import Namespace
import os
from re import findall, MULTILINE
from StringIO import StringIO
import sys
from testfixtures import TempDirectory
import unittest
import glob

import jacquard.utils as utils
import jacquard.tag as tag
from jacquard.vcf import FileReader
import jacquard.logger as logger
import jacquard.variant_callers.strelka as strelka
import jacquard.variant_callers.varscan as varscan

mock_log_called = False
mock_log_messages = []

class MockVcfRecord(object):
    def __init__(self, ref = "foo", alt = "bar", filter = "baz"):
        self.ref = "foo"
        self.alt = "bar"
        self.filter = "baz"
        self.content = "foo"
  
mock_vcf_record = MockVcfRecord()

def mock_log(msg, *args):
    global mock_log_called
    mock_log_called = True
    global mock_log_messages
    mock_log_messages.append(msg.format(*[str(i) for i in args]))

class MockWriter():
    def __init__(self):
        self._content = []
        self.output_filepath = "foo"
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
        return vcfRecord.content

    def get_new_metaheaders(self):
        return self.metaheaders

class MockVcfReader(object):
    def __init__(self, input_filepath="vcfName", metaheaders=["##metaheaders"], column_header="#header", file_name = "foo"):
        self.input_filepath = input_filepath
        self.file_name = file_name
        self.metaheaders = metaheaders
        self.caller = MockCaller()
        self.column_header = column_header
        self.opened = False
        self.closed = False

    def open(self):
        self.opened = True

    def vcf_records(self):
        global mock_vcf_record
        yield mock_vcf_record

    def close(self):
        self.closed = True

def build_mock_get_caller_method(callers):
    def get_caller(metaheaders, column_header, name):
        if len(callers)>1:
            for caller in callers:
                if caller.name == name:
                    return caller
        else:
            return callers[0]
    return get_caller

class TagTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output
        self.original_info = logger.info
        self.original_error = logger.error
        self.original_warning = logger.warning
        self.original_debug = logger.debug
        self._change_mock_logger()

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr
        self._reset_mock_logger()

    def _change_mock_logger(self):
        global mock_log_called
        mock_log_called = False
        global mock_log
        logger.info = mock_log
        logger.error = mock_log
        logger.warning = mock_log
        logger.debug = mock_log

    def _reset_mock_logger(self):
        logger.info = self.original_info
        logger.error = self.original_error
        logger.warning = self.original_warning
        logger.debug = self.original_debug
        global mock_log_messages
        mock_log_messages = []

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

            self.assertEqual("A.vcf", vcf_readers[0].file_name)
            self.assertEqual(["##source=strelka"], vcf_readers[0].metaheaders)
            self.assertEqual("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR",
                              vcf_readers[0].column_header)
            self.assertEqual("B.vcf", vcf_readers[1].file_name)
            self.assertEqual(["##source=strelka"], vcf_readers[1].metaheaders)
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

            actual_readers = sorted([reader.file_name for reader in actual_dict])
            expected_readers = ["A.vcf","B.vcf"]
            self.assertEquals(actual_readers, expected_readers)

            actual_writers = sorted([writer.output_filepath for writer in actual_dict.values()])
            expected_writers = [os.path.join(output_dir.path, base_filename) for base_filename in ["A.jacquardTags.vcf", "B.jacquardTags.vcf"]]
            self.assertEquals(actual_writers, expected_writers)

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
            self.assertEquals(1, len(output_lines))
            global mock_log_called
            self.assertTrue(mock_log_called)
#             self.assertRegexpMatches(output_lines[0], 'Recognized \[2\] Strelka file\(s\)')

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

            self.assertRaisesRegexp(utils.JQException,
                                    "VCF files could not be parsed",
                                    tag._build_vcf_readers,
                                    input_dir.path)

    def test_tag_files(self):
        reader = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader")
        reader.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        writer = MockWriter()
        vcf_readers_to_writers = {reader: writer}
        execution_context = []
        tag.tag_files(vcf_readers_to_writers, execution_context, build_mock_get_caller_method([MockCaller()]))

        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertEquals(["##originalMeta1",
                           "##originalMeta2",
                           "##jacquard.tag.caller=MockCaller",
                           "##mockCallerMetaheader1",
                           '##FILTER=<ID=JQ_EXCLUDE,Description="This variant record is problematic and will be excluded from downstream Jacquard processing.",Source="Jacquard",Version="">',
                           '##FILTER=<ID=JQ_MALFORMED_REF,Description="The format of the reference value for this variant record does not comply with VCF standard.",Source="Jacquard",Version="">',
                           '##FILTER=<ID=JQ_MALFORMED_ALT,Description="The the format of the alternate allele value for this variant record does not comply with VCF standard.",Source="Jacquard",Version="">',
                           "#columnHeader",
                           "foo"], writer.lines())

        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)

    def test_tag_files_proper(self):
        class MockVcfRecord(object):
            def __init__(self):
                self.ref = "A"
                self.alt = "A"
                self.filter = "filter"
                self.content = "foo"

        global mock_vcf_record
        mock_vcf_record = MockVcfRecord()

        reader = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader")
        reader.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        writer = MockWriter()

        vcf_readers_to_writers = {reader: writer}
        execution_context = []
        tag.tag_files(vcf_readers_to_writers, execution_context, build_mock_get_caller_method([MockCaller()]))

        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertEquals(["##originalMeta1",
                           "##originalMeta2",
                           "##jacquard.tag.caller=MockCaller",
                           "##mockCallerMetaheader1",
                           "#columnHeader",
                           "foo"], writer.lines())

        self.assertEquals("filter", mock_vcf_record.filter)
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)

    def test_tag_files_proper_multipleFilters(self):
        class MockVcfRecord(object):
            def __init__(self):
                self.ref = "A"
                self.alt = "A"
                self.filter = "filter;QSS_ref"
                self.content = "foo"

        global mock_vcf_record
        mock_vcf_record = MockVcfRecord()

        reader = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader")
        reader.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        writer = MockWriter()

        vcf_readers_to_writers = {reader: writer}
        execution_context = []
        tag.tag_files(vcf_readers_to_writers, execution_context, build_mock_get_caller_method([MockCaller()]))

        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertEquals(["##originalMeta1",
                           "##originalMeta2",
                           "##jacquard.tag.caller=MockCaller",
                           "##mockCallerMetaheader1",
                           "#columnHeader",
                           "foo"], writer.lines())

        self.assertEquals("filter;QSS_ref", mock_vcf_record.filter)
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)

    def test_tag_files_malformedRef(self):
        class MockVcfRecord(object):
            def __init__(self):
                self.ref = "/"
                self.alt = "A"
                self.filter = "filter"
                self.content = "foo"

        global mock_vcf_record
        mock_vcf_record = MockVcfRecord()

        reader = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader")
        reader.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        writer = MockWriter()

        vcf_readers_to_writers = {reader: writer}
        execution_context = []
        tag.tag_files(vcf_readers_to_writers, execution_context, build_mock_get_caller_method([MockCaller()]))

        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertEquals(["##originalMeta1",
                           "##originalMeta2",
                           "##jacquard.tag.caller=MockCaller",
                           "##mockCallerMetaheader1",
                           '##FILTER=<ID=JQ_EXCLUDE,Description="This variant record is problematic and will be excluded from downstream Jacquard processing.",Source="Jacquard",Version="">',
                           '##FILTER=<ID=JQ_MALFORMED_REF,Description="The format of the reference value for this variant record does not comply with VCF standard.",Source="Jacquard",Version="">',
                           "#columnHeader",
                           "foo"], writer.lines())

        self.assertEquals("filter;JQ_EXCLUDE;JQ_MALFORMED_REF", mock_vcf_record.filter)
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)

        self.assertIn("foo|Added filter flag [JQ_MALFORMED_REF] to [1] variant records.", mock_log_messages)
        self.assertIn("Added a filter flag to [1] problematic MockCaller variant records.", mock_log_messages)
        self.assertIn("A total of [1] problematic variant records failed Jacquard's filters. See output and log for details.", mock_log_messages)

    def test_tag_files_multipleCallers(self):
        class MockVcfRecord(object):
            def __init__(self):
                self.ref = "/"
                self.alt = "A"
                self.filter = "filter"
                self.content = "foo"

        global mock_vcf_record
        mock_vcf_record = MockVcfRecord()

        reader1 = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader", file_name="foo")
        reader2 = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader", file_name="bar")
        reader1.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        reader2.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        writer1 = MockWriter()
        writer2 = MockWriter()

        vcf_readers_to_writers = {reader1: writer1, reader2: writer2}
        execution_context = []
        tag.tag_files(vcf_readers_to_writers, execution_context,
                      build_mock_get_caller_method([MockCaller(name="foo"),
                                                    MockCaller(name="bar")]))

        self.assertIn("foo|Added filter flag [JQ_MALFORMED_REF] to [1] variant records.", mock_log_messages)
        self.assertIn("Added a filter flag to [1] problematic foo variant records.", mock_log_messages)
        self.assertIn("Added a filter flag to [1] problematic bar variant records.", mock_log_messages)
        self.assertIn("A total of [2] problematic variant records failed Jacquard's filters. See output and log for details.", mock_log_messages)

    def test_tag_files_malformedRef_emptyFilter(self):
        class MockVcfRecord(object):
            def __init__(self):
                self.ref = "/"
                self.alt = "A"
                self.filter = "."
                self.content = "foo"

        global mock_vcf_record
        mock_vcf_record = MockVcfRecord()

        reader = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader")
        reader.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        writer = MockWriter()

        vcf_readers_to_writers = {reader: writer}
        execution_context = []
        tag.tag_files(vcf_readers_to_writers, execution_context, build_mock_get_caller_method([MockCaller()]))

        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertEquals(["##originalMeta1",
                           "##originalMeta2",
                           "##jacquard.tag.caller=MockCaller",
                           "##mockCallerMetaheader1",
                           '##FILTER=<ID=JQ_EXCLUDE,Description="This variant record is problematic and will be excluded from downstream Jacquard processing.",Source="Jacquard",Version="">',
                           '##FILTER=<ID=JQ_MALFORMED_REF,Description="The format of the reference value for this variant record does not comply with VCF standard.",Source="Jacquard",Version="">',
                           "#columnHeader",
                           "foo"], writer.lines())

        self.assertEquals("JQ_EXCLUDE;JQ_MALFORMED_REF",mock_vcf_record.filter)
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)

        self.assertIn("foo|Added filter flag [JQ_MALFORMED_REF] to [1] variant records.", mock_log_messages)
        self.assertIn("Added a filter flag to [1] problematic MockCaller variant records.", mock_log_messages)
        self.assertIn("A total of [1] problematic variant records failed Jacquard's filters. See output and log for details.", mock_log_messages)

    def test_tag_files_malformedAlt(self):
        class MockVcfRecord(object):
            def __init__(self):
                self.ref = "A"
                self.alt = "/"
                self.filter = "filter"
                self.content = "foo"

        global mock_vcf_record
        mock_vcf_record = MockVcfRecord()

        reader = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader")
        reader.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        writer = MockWriter()

        vcf_readers_to_writers = {reader: writer}
        execution_context = []
        tag.tag_files(vcf_readers_to_writers, execution_context, build_mock_get_caller_method([MockCaller()]))

        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertEquals(["##originalMeta1",
                           "##originalMeta2",
                           "##jacquard.tag.caller=MockCaller",
                           "##mockCallerMetaheader1",
                           '##FILTER=<ID=JQ_EXCLUDE,Description="This variant record is problematic and will be excluded from downstream Jacquard processing.",Source="Jacquard",Version="">',
                           '##FILTER=<ID=JQ_MALFORMED_ALT,Description="The the format of the alternate allele value for this variant record does not comply with VCF standard.",Source="Jacquard",Version="">',
                           "#columnHeader",
                           "foo"], writer.lines())

        self.assertEquals("filter;JQ_EXCLUDE;JQ_MALFORMED_ALT",mock_vcf_record.filter)
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)

        self.assertIn("foo|Added filter flag [JQ_MALFORMED_ALT] to [1] variant records.", mock_log_messages)
        self.assertIn("Added a filter flag to [1] problematic MockCaller variant records.", mock_log_messages)
        self.assertIn("A total of [1] problematic variant records failed Jacquard's filters. See output and log for details.", mock_log_messages)

    def test_tag_files_missingAltDot(self):
        class MockVcfRecord(object):
            def __init__(self):
                self.ref = "A"
                self.alt = "."
                self.filter = "filter"
                self.content = "foo"

        global mock_vcf_record
        mock_vcf_record = MockVcfRecord()

        reader = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader")
        reader.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        writer = MockWriter()

        vcf_readers_to_writers = {reader: writer}
        execution_context = []
        tag.tag_files(vcf_readers_to_writers, execution_context, build_mock_get_caller_method([MockCaller()]))

        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertEquals(["##originalMeta1",
                           "##originalMeta2",
                           "##jacquard.tag.caller=MockCaller",
                           "##mockCallerMetaheader1",
                           '##FILTER=<ID=JQ_EXCLUDE,Description="This variant record is problematic and will be excluded from downstream Jacquard processing.",Source="Jacquard",Version="">',
                           '##FILTER=<ID=JQ_MISSING_ALT,Description="The alternate allele is missing for this variant record.",Source="Jacquard",Version="">',
                           "#columnHeader",
                           "foo"], writer.lines())

        self.assertEquals("filter;JQ_EXCLUDE;JQ_MISSING_ALT",mock_vcf_record.filter)
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)

        self.assertIn("foo|Added filter flag [JQ_MISSING_ALT] to [1] variant records.", mock_log_messages)
        self.assertIn("Added a filter flag to [1] problematic MockCaller variant records.", mock_log_messages)
        self.assertIn("A total of [1] problematic variant records failed Jacquard's filters. See output and log for details.", mock_log_messages)

    def test_tag_files_missingAltAsterisk(self):
        class MockVcfRecord(object):
            def __init__(self):
                self.ref = "A"
                self.alt = "*"
                self.filter = "filter"
                self.content = "foo"

        global mock_vcf_record
        mock_vcf_record = MockVcfRecord()

        reader = MockVcfReader(metaheaders=["##originalMeta1", "##originalMeta2"], column_header="#columnHeader")
        reader.caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        writer = MockWriter()

        vcf_readers_to_writers = {reader: writer}
        execution_context = []
        tag.tag_files(vcf_readers_to_writers, execution_context, build_mock_get_caller_method([MockCaller()]))

        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertEquals(["##originalMeta1",
                           "##originalMeta2",
                           "##jacquard.tag.caller=MockCaller",
                           "##mockCallerMetaheader1",
                           '##FILTER=<ID=JQ_EXCLUDE,Description="This variant record is problematic and will be excluded from downstream Jacquard processing.",Source="Jacquard",Version="">',
                           '##FILTER=<ID=JQ_MISSING_ALT,Description="The alternate allele is missing for this variant record.",Source="Jacquard",Version="">',
                           "#columnHeader",
                           "foo"], writer.lines())

        self.assertEquals("filter;JQ_EXCLUDE;JQ_MISSING_ALT",mock_vcf_record.filter)
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)

        self.assertIn("foo|Added filter flag [JQ_MISSING_ALT] to [1] variant records.", mock_log_messages)
        self.assertIn("Added a filter flag to [1] problematic MockCaller variant records.", mock_log_messages)
        self.assertIn("A total of [1] problematic variant records failed Jacquard's filters. See output and log for details.", mock_log_messages)

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
        vcf_content1 = vcf_content1.replace('|', "\t")
        vcf_content2 = vcf_content2.replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.vcf", vcf_content1)
            input_dir.write("B.vcf", vcf_content2)

            args = Namespace(input=input_dir.path,
                             output=output_dir.path)

            tag.execute(args, [])

            output_dir.check("A.jacquardTags.vcf", "B.jacquardTags.vcf")
            file_content1 = output_dir.read("A.jacquardTags.vcf")
            file_content2 = output_dir.read("B.jacquardTags.vcf")

            self.assertTrue('##FORMAT=<ID={0}HC_SOM'.format(strelka.JQ_STRELKA_TAG) in file_content1)
            self.assertEquals(2,
                              len(findall(r'^chr.*{0}HC_SOM'.format(strelka.JQ_STRELKA_TAG), file_content1, MULTILINE)))
            self.assertTrue('##FORMAT=<ID={0}HC_SOM'.format(varscan.JQ_VARSCAN_TAG) in file_content2)
            self.assertEquals(2,
                              len(findall(r'^chr.*{0}HC_SOM'.format(varscan.JQ_VARSCAN_TAG), file_content2, MULTILINE)))

    def test_functional_tag(self):
        with TempDirectory() as output_dir:
            module_testdir = os.path.dirname(os.path.realpath(__file__))+"/functional_tests/02_tag"
            input_dir = os.path.join(module_testdir,"input")
            args = Namespace(input=input_dir,
                         output=output_dir.path)

            execution_context = ["##jacquard.version={0}".format(utils.__version__),
                "##jacquard.command={0}".format(" ".join(["tag",
                                        os.path.join(module_testdir,"input"),
                                        os.path.join(module_testdir,"benchmark")])),
                "##jacquard.cwd="
                ]
            tag.execute(args,execution_context)

            output_file = glob.glob(os.path.join(output_dir.path, "*.vcf"))[0]

            actual_file = FileReader(os.path.join(output_dir.path, output_file))
            actual_file.open()
            actual = []
            for line in actual_file.read_lines():
                actual.append(line)
            actual_file.close()

            module_outdir = os.path.join(module_testdir, "benchmark")
            output_file = os.listdir(module_outdir)[0]
            expected_file = FileReader(os.path.join(module_outdir, output_file))
            expected_file.open()
            expected = []
            for line in expected_file.read_lines():
                expected.append(line)
            expected_file.close()

            self.assertEquals(len(expected), len(actual))
            self.assertEquals(133, len(actual))

            for i in xrange(len(expected)):
                if expected[i].startswith("##jacquard.cwd="):
                    self.assertTrue(actual[i].rstrip() in expected[i])
                elif expected[i].startswith("##jacquard.command="):
                    self.assertTrue(actual[i].startswith("##jacquard.command="))
                else:
                    self.assertEquals(expected[i].rstrip(), actual[i].rstrip())

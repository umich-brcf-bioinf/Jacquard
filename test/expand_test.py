# pylint: disable=C0103,C0301,R0903,R0904,W0603,W0613,W0212,C0111
from argparse import Namespace
from collections import OrderedDict
import glob
import os
from testfixtures import TempDirectory
import unittest

import jacquard.utils as utils
import jacquard.logger as logger
import jacquard.expand as expand
from jacquard.vcf import FileReader, VcfReader

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
mock_log_called = False
mock_warnings = []

def mock_log(msg, *args):
    global mock_log_called
    mock_log_called = True

def mock_warning(msg, *args):
    global mock_warnings
    mock_warnings.append(msg.format(*[str(i) for i in args]))

def _change_mock_logger():
    global mock_log_called
    mock_log_called = False
    global mock_warnings
    mock_warnings = []
    logger.info = mock_log
    logger.error = mock_log
    logger.warning = mock_warning
    logger.debug = mock_log

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

# pylint: disable=W0102
class MockVcfReader(object):
    def __init__(self,
                 input_filepath="vcfName",
                 metaheaders=["##metaheaders"],
                 column_header="#header",
                 content=["foo"]):
        self.input_filepath = input_filepath
        self.metaheaders = metaheaders
        self.column_header = column_header
        self.opened = False
        self.closed = False
        self.content = content

    def open(self):
        self.opened = True

    def vcf_records(self):
        for content in self.content:
            yield MockVcfRecord(content)

    def close(self):
        self.closed = True

class MockVcfRecord(object):
    def __init__(self, content):
        self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, \
            self.filter, self.info, self.format = content[0:9]
        self.samples = content[9:]

        tags = self.format.split(":")
        self.sample_dict = {}
        for i, sample in enumerate(self.samples):
            values = sample.split(":")
            self.sample_dict[i] = OrderedDict(zip(tags, values))

    def get_info_dict(self):
        info_dict = {}

        for key_value in self.info.split(";"):
            if "=" in key_value:
                key, value = key_value.split("=")
                info_dict[key] = value
            else:
                info_dict[key_value] = key_value

        return info_dict

class MockFileWriter(object):
    def __init__(self):
        self.written = []

    def write(self, text):
        self.written.append(text)

class ExpandTestCase(unittest.TestCase):
    def setUp(self):
        self.original_info = logger.info
        self.original_error = logger.error
        self.original_warning = logger.warning
        self.original_debug = logger.debug
        _change_mock_logger()

    def tearDown(self):
        self._reset_mock_logger()

    def _reset_mock_logger(self):
        logger.info = self.original_info
        logger.error = self.original_error
        logger.warning = self.original_warning
        logger.debug = self.original_debug

    def test_parse_meta_headers(self):
        meta_headers = ['##ALT=<ID=DEL,Description="Deletion">',
                        '##INFO=<ID=AC,Number=.,Description="foo">]',
                        '##INFO=<ID=AA,Number=1,Description="Ancestral Allele"',
                        '##FORMAT=<ID=SP,Type=Integer,Description="bar">',
                        '##RUNTIME_ARG=allele freq. cutoff: 5']
        (info_fields, format_tags) = expand._parse_meta_headers(meta_headers)

        self.assertEquals(["AA", "AC"], info_fields)
        self.assertEquals(["SP"], format_tags)

    def test_parse_meta_headers_missing(self):
        meta_headers = ['##ALT=<ID=DEL,Description="Deletion">',
                        '##INFO=<ID=AC,Number=.,Description="foo">]',
                        '##INFO=<ID=AA,Number=1,Description="Ancestral Allele"',
                        '##RUNTIME_ARG=allele freq. cutoff: 5']

        self.assertRaisesRegexp(utils.JQException,
                                ("Unable to parse meta_headers for INFO "
                                 "and/or FORMAT fields. Review input and "
                                 "try again."),
                                expand._parse_meta_headers,
                                meta_headers)

    def test_create_row_dict(self):
        column_list = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                       "INFO", "FORMAT", "SAMPLE_A|NORMAL", "SAMPLE_A|TUMOR"]
        row = ["1", "42", "rs32", "A", "AT", "30", "PASS",
               "SNP;SOMATIC=1", "DP:AF", "50:0.2", "87:0.3"]
        vcf_record = MockVcfRecord(row)

        actual_dict = expand._create_row_dict(column_list, vcf_record)

        expected_dict = {"CHROM": "1",
                         "POS": "42",
                         "ID": "rs32",
                         "REF": "A",
                         "ALT": "AT",
                         "QUAL": "30",
                         "FILTER": "PASS",
                         "SNP": "SNP",
                         "SOMATIC": "1",
                         "DP|SAMPLE_A|NORMAL": "50",
                         "DP|SAMPLE_A|TUMOR": "87",
                         "AF|SAMPLE_A|NORMAL": "0.2",
                         "AF|SAMPLE_A|TUMOR": "0.3"}
        self.assertEquals(expected_dict, actual_dict)

    def test_create_actual_column_list(self):
        potential_col_list = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                              "FILTER", "MULT_ALT", "DP|SAMPLE_A|NORMAL",
                              "AF|SAMPLE_A|NORMAL", "DP|SAMPLE_B|NORMAL",
                              "AF|SAMPLE_B|NORMAL"]
        col_spec = ["CHROM", r"DP\|.*", "POS"]
        actual_column_list = expand._create_actual_column_list(col_spec,
                                                               potential_col_list,
                                                               "col_spec.txt")
        expected_column_list = ["CHROM",
                                "DP|SAMPLE_A|NORMAL",
                                "DP|SAMPLE_B|NORMAL",
                                "POS"]

        self.assertEquals(expected_column_list, actual_column_list)

    def test_create_actual_column_list_noColsIncluded(self):
        potential_col_list = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL"]
        col_spec = ["FOO", "BAR", "BAZ"]

        self.assertRaises(utils.JQException, 
                          expand._create_actual_column_list,
                          col_spec, 
                          potential_col_list,
                          "col_spec.txt")

    def test_create_actual_column_list_regexNotInCols(self):
        potential_col_list = ["CHROM", "POS", "ID", "REF", "ALT"]
        col_spec = ["CHROM", "FOO", "POS", "BAR"]
        col_spec_fname = "spec.txt"
        actual_column_list = expand._create_actual_column_list(col_spec,
                                                               potential_col_list,
                                                               col_spec_fname)
        self.assertEquals(["CHROM", "POS"], actual_column_list)

        exp_warning_1 = ("The expression [FOO] in column specification "+
                         "file [spec.txt:2] didn't match any input columns; "+
                         "columns may have matched earlier expressions, or "+
                         "this expression may be irrelevant.")
        exp_warning_2 = ("The expression [BAR] in column specification "+
                         "file [spec.txt:4] didn't match any input columns; "+
                         "columns may have matched earlier expressions, or "+
                         "this expression may be irrelevant.")

        self.assertEquals([exp_warning_1, exp_warning_2], mock_warnings)

    def test_create_potential_column_list(self):
        file_contents = ['##INFO=<ID=AF,Number=1>\n',
                         '##INFO=<ID=AA,Number=1>\n',
                         '##FORMAT=<ID=GT,Number=1>\n',
                         '##FORMAT=<ID=GQ,Number=1,Description="bar">\n',
                         '#chrom\tpos\tid\tref\talt\n',
                         'record1\n',
                         'record2']
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        vcf_reader = VcfReader(mock_file_reader)

        actual_col_list = expand._create_potential_column_list(vcf_reader)
        expected_col_list = ["chrom", "pos", "id", "ref", "alt",
                             "AF", "AA",
                             "GT", "GQ"]
        self.assertEquals(expected_col_list, actual_col_list)

    def test_disambiguate_column_names(self):
        column_header = ["CHROM", "POS", "ID", "REF"]
        info_header = ["HOM", "AA", "SOM"]

        actual = expand._disambiguate_column_names(column_header, info_header)
        expected = ["HOM", "AA", "SOM"]

        self.assertEquals(expected, actual)

        column_header = ["CHROM", "POS", "ID", "REF"]
        info_header = ["HOM", "AA", "ID", "SOM"]
        actual = expand._disambiguate_column_names(column_header, info_header)
        expected = ["INFO_HOM", "INFO_AA", "INFO_ID", "INFO_SOM"]

        self.assertEquals(expected, actual)

    def test_append_format_tags_to_samples(self):
        format_tags = ["foo", "bar"]
        samples = ["sampleB", "sampleA"]
        actual = expand._append_format_tags_to_samples(format_tags, samples)

        expected = ["bar|sampleA", "bar|sampleB", "foo|sampleA", "foo|sampleB"]

        self.assertEquals(expected, actual)

    def test_get_headers(self):
        meta_headers = ['##ALT=<ID=DEL,Description="Deletion">',
                        '##INFO=<ID=AC,Number=.,Description="foo">]',
                        '##INFO=<ID=AA,Number=1,Description="Ancestral Allele"',
                        '##FORMAT=<ID=SP,Number=1,Description="bar">',
                        '##RUNTIME_ARG=allele freq. cutoff: 5']

        col_header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\tsampleB"
        mock_reader = MockVcfReader(metaheaders=meta_headers,
                                    column_header=col_header)
        actual = expand._get_headers(mock_reader)

        expected = (["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"],
                    ["AA", "AC"],
                    ["SP|sampleA", "SP|sampleB"])

        self.assertEquals(expected, actual)

    def test_write_vcf_records(self):
        column_header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA"

        mock_vcf_reader = MockVcfReader(content=[["CHROM",
                                                  "POS",
                                                  "ID",
                                                  "REF",
                                                  "ALT",
                                                  "QUAL",
                                                  "FILTER",
                                                  "tag1=val1;tag3=val3;tag4",
                                                  "FOO:BAR",
                                                  "42:1"]],
                                        column_header=column_header)

        mock_file_writer = MockFileWriter()

        info_header = ["tag1", "tag2", "tag3", "tag4"]
        format_sample_header = ["BAR|sampleA", "FOO|sampleA"]

        split_column_header = column_header.split("\t")[0:7]

        header_dict = OrderedDict([("column_header", split_column_header),
                                   ("info_header", info_header),
                                   ("format_header", format_sample_header)])

        expand._write_vcf_records(mock_vcf_reader, mock_file_writer, header_dict)

        actual = mock_file_writer.written
        expected = ["CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tval1\t\tval3\ttag4\t1\t42\n"]

        self.assertEquals(expected, actual)

    def test_filter_and_sort(self):
        header = OrderedDict([("column_header", ["CHROM", "POS", "ID"]),
                              ("info_header", ["infoA", "infoB", "infoC"]),
                              ("format_header", ["tagA|sample1",
                                                 "tagB|sample1",
                                                 "tagC|sample1"])])
        columns_to_expand = ["CHROM", "infoA", r"tagA\|.*", "infoB"]
        actual_header_dict = expand._filter_and_sort(header, columns_to_expand)
        expected_header = {'column_header': ["CHROM"],
                           'format_header': ['tagA|sample1'],
                           'info_header': ["infoA", "infoB"]}

        self.assertEquals(expected_header, actual_header_dict)

    def test_filter_and_sort_matchesBeginningOfColumnByDefault(self):
        header = OrderedDict([("column_header", ["POS"]),
                              ("info_header", ["A_POS", "A_POS_B", "POS_B"])])

        actual_header_dict = expand._filter_and_sort(header, [r"POS"])
        expected_header = {'column_header': ["POS"]}
        self.assertEquals(actual_header_dict, expected_header)

        actual_header_dict = expand._filter_and_sort(header, [r".*POS.*"])
        expected_header = {'column_header': ["POS"],
                           'info_header': ["A_POS", "A_POS_B", "POS_B"]}
        self.assertEquals(actual_header_dict, expected_header)

        actual_header_dict = expand._filter_and_sort(header, [r"^POS.*"])
        expected_header = {'column_header': ["POS"],
                           'info_header': ["POS_B"]}
        self.assertEquals(actual_header_dict, expected_header)

    def test_filter_and_sort_missingInfo(self):
        header = OrderedDict([("column_header", ["CHROM", "POS", "ID"]),
                              ("info_header", ["infoA", "infoB", "infoC"]),
                              ("format_header", ["tagA|sample1",
                                                 "tagB|sample1",
                                                 "tagC|sample1"])])

        columns_to_expand = [r"CHROM", r"tagA\|.*"]
        actual_header = expand._filter_and_sort(header, columns_to_expand)

        expected_header = {'column_header': ["CHROM"],
                           'format_header': ['tagA|sample1']}

        self.assertEquals(expected_header, actual_header)

    def test_filter_and_sort_checkOrder(self):
        header = OrderedDict([("column_header", ["CHROM", "POS", "ID", "REF"]),
                              ("info_header", ["infoA", "infoB", "infoC"]),
                              ("format_header", ["tagA|sample1",
                                                 "tagB|sample1",
                                                 "tagC|sample1"])])

        columns_to_expand = [r"CHROM", r"tagA\|.*", r"ID", r"POS"]
        actual_header = expand._filter_and_sort(header, columns_to_expand)

        expected_header = {'column_header': ["CHROM", "ID", "POS"],
                           'format_header': ['tagA|sample1']}

        self.assertEquals(expected_header, actual_header)

    def test_filter_and_sort_noColumsIncluded(self):
        header = OrderedDict([("column_header", ["CHROM", "POS", "ID", "REF"]),
                              ("info_header", ["infoA", "infoB", "infoC"]),
                              ("format_header", ["tagA|sample1",
                                                 "tagB|sample1",
                                                 "tagC|sample1"])])

        columns_to_expand = ["^foo$", "^bar*"]

        self.assertRaisesRegexp(utils.JQException,
                                ("The column specification file would "
                                 "exclude all input columns. Review inputs/"
                                 "usage and try again."),
                                expand._filter_and_sort, header, columns_to_expand)

    def test_execute_files(self):
        vcf_content = ('''##source=strelka
##FORMAT=foo
##INFO=bar
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
''').replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("P1.vcf", vcf_content)
            input_file = os.path.join(input_dir.path, "P1.vcf")
            output_file = os.path.join(output_dir.path, "P1.txt")
            args = Namespace(input=input_file,
                             output=output_file,
                             column_specification=0)

            expand.execute(args, ["extra_header1", "extra_header2"])

            output_dir.check("P1.txt")
            with open(os.path.join(output_dir.path, "P1.txt")) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()

        self.assertEquals(3, len(actual_output_lines))

    def test_execute_dirs(self):
        vcf_content = ('''##source=strelka
##FORMAT=foo
##INFO=bar
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
''').replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("P1.vcf", vcf_content)
            input_dir.write("P2.vcf", vcf_content)
            args = Namespace(input=input_dir.path,
                             output=output_dir.path,
                             column_specification=0)

            expand.execute(args, ["extra_header1", "extra_header2"])

            output_dir.check("P1.txt", "P2.txt")
            with open(os.path.join(output_dir.path, "P1.txt")) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()

        self.assertEquals(3, len(actual_output_lines))

    def test_expand_emptyInputDir(self):

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            args = Namespace(input=input_dir.path,
                             output=output_dir.path,
                             column_specification=0)

            self.assertRaisesRegexp(utils.JQException,
                                    ("Specified input directory .* contains "
                                     "no VCF files. Review inputs and .*"),
                                    expand.execute,
                                    args,
                                    ["extra_header1", "extra_header2"])

    def test_execute_inputFileOutputDir(self):
        vcf_content = ('''##source=strelka
##FORMAT=foo
##INFO=bar
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
''').replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("P1.vcf", vcf_content)
            input_file = os.path.join(input_dir.path, "P1.vcf")
            args = Namespace(input=input_file,
                             output=output_dir.path,
                             column_specification=0)

            self.assertRaisesRegexp(utils.JQException,
                                    ("Specified output .* must be a file if "
                                     "input .* is a file."),
                                    expand.execute,
                                    args,
                                    ["extra_header1", "extra_header2"])

    def test_execute_inputDirOutputFile(self):
        vcf_content = ('''##source=strelka
##FORMAT=foo
##INFO=bar
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR
''').replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("P1.vcf", vcf_content)
            output_dir.write("P2.vcf", vcf_content)
            output_file = os.path.join(output_dir.path, "P2.vcf")
            args = Namespace(input=input_dir.path,
                             output=output_file,
                             column_specification=0)

            self.assertRaisesRegexp(utils.JQException,
                                    ("Specified output .* must be a directory "
                                     "if input .* is a directory."),
                                    expand.execute,
                                    args,
                                    ["extra_header1", "extra_header2"])

    def test_execute_colSpecValid(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir, TempDirectory() as col_spec_dir:
            col_spec_dir.write("col_spec.txt", "foo\nbar")
            col_spec_file = os.path.join(col_spec_dir.path, "col_spec.txt")
            args = Namespace(input=input_dir.path,
                             output=output_dir.path,
                             column_specification=col_spec_file)

            self.assertRaisesRegexp(utils.JQException,
                                    ("Specified input directory .* contains "
                                     "no VCF files. Review inputs .*"),
                                    expand.execute,
                                    args,
                                    ["extra_header1", "extra_header2"])

    def test_execute_colSpecInvalid(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir, TempDirectory() as col_spec_dir:
            col_spec_dir.write("col_spec.txt", "foo\nbar")
            args = Namespace(input=input_dir.path,
                             output=output_dir.path,
                             column_specification=col_spec_dir.path)

            self.assertRaisesRegexp(utils.JQException, "The column "
                                    "specification file .* could not be "
                                    "read. Review inputs/usage and try again",
                                    expand.execute, args,
                                    ["extra_header1", "extra_header2"])

    def test_functional_expand(self):
        with TempDirectory() as output_dir:
            module_testdir = os.path.dirname(os.path.realpath(__file__))+"/functional_tests/06_expand"
            input_dir = os.path.join(module_testdir,"input")
            args = Namespace(input=input_dir,
                         output=output_dir.path,
                         column_specification=None)

            execution_context = ["##jacquard.version={0}".format(utils.__version__),
                "##jacquard.command=",
                "##jacquard.cwd="]
            expand.execute(args,execution_context)

            output_file = glob.glob(os.path.join(output_dir.path, "consensus.txt"))[0]

            actual_file = FileReader(output_file)
            actual_file.open()
            actual = []
            for line in actual_file.read_lines():
                actual.append(line)
            actual_file.close()

            module_outdir = os.path.join(module_testdir,"benchmark")
            output_file = os.listdir(module_outdir)[0]
            expected_file = FileReader(os.path.join(module_outdir,output_file))
            expected_file.open()
            expected = []
            for line in expected_file.read_lines():
                expected.append(line)
            expected_file.close()

            self.assertEquals(len(expected), len(actual))

            self.assertEquals(11, len(actual))

            for i in xrange(len(expected)):
                if expected[i].startswith("##jacquard.cwd="):
                    self.assertTrue(actual[i].startswith("##jacquard.cwd="))
                elif expected[i].startswith("##jacquard.command="):
                    self.assertTrue(actual[i].startswith("##jacquard.command="))
                else:
                    self.assertEquals(expected[i], actual[i])


#pylint: disable=line-too-long, too-many-instance-attributes, unused-argument
#pylint: disable=global-statement, star-args, too-few-public-methods, invalid-name
#pylint: disable=too-many-public-methods
from __future__ import absolute_import
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
import test.test_case as test_case

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
MOCK_LOG_CALLED = False
MOCK_WARNINGS = []

def mock_log(msg, *args):
    global MOCK_LOG_CALLED
    MOCK_LOG_CALLED = True

def mock_warning(msg, *args):
    MOCK_WARNINGS.append(msg.format(*[str(i) for i in args]))

def _change_mock_logger():
    global MOCK_LOG_CALLED
    MOCK_LOG_CALLED = False
    global MOCK_WARNINGS
    MOCK_WARNINGS = []
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

        self.assertEquals([exp_warning_1, exp_warning_2], MOCK_WARNINGS)

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
                             "AA", "AF",
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

    def test_execute_files(self):
        vcf_content = ('''##source=strelka
##INFO=<ID=SOMATIC,Number=1,Description="foo">
##FORMAT=<ID=GT,Number=1,Description="bar">
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|SOMATIC|GT|0/1|0/1
chr2|1|.|A|C|.|.|SOMATIC|GT|0/1|0/1
''').replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("P1.vcf", vcf_content)
            input_file = os.path.join(input_dir.path, "P1.vcf")
            output_file = os.path.join(output_dir.path, "P1.txt")
            args = Namespace(input=input_file,
                             output=output_file,
                             column_specification=0)

            expand.execute(args, ["##extra_header1", "##extra_header2"])

            output_dir.check("P1.txt")
            with open(os.path.join(output_dir.path, "P1.txt")) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()
        print actual_output_lines
        self.assertEquals(3, len(actual_output_lines))

    def test_execute_dirs(self):
        vcf_content = ('''##source=strelka
##INFO=<ID=SOMATIC,Number=1,Description="foo">
##FORMAT=<ID=GT,Number=1,Description="bar">
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

            expand.execute(args, ["##extra_header1", "##extra_header2"])

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
##INFO=<ID=SOMATIC,Number=1,Description="foo">
##FORMAT=<ID=GT,Number=1,Description="bar">
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
##INFO=<ID=AA,Number=1,Description="Ancestral Allele"'
##FORMAT=<ID=SP,Number=1,Description="bar">'r
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

##TODO: fix this test!
    def xtest_functional_expand(self):
        #these meta_headers break it (among others) because of commas:
        ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele, ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/technical/reference/ancestral_alignments/README">
        ##INFO=<ID=COSMIC,Number=.,Type=String,Description="Data from COSMIC. Format: COSMIC=(seqname|source|feature|start|end|score|strand|frame|attributes|comments),",Source=Epee|COSMIC>
        ##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">

        with TempDirectory() as output_dir:
            module_testdir = os.path.dirname(os.path.realpath(__file__))+"/functional_tests/06_expand"
            input_dir = os.path.join(module_testdir,"input")
            args = Namespace(input=input_dir,
                         output=output_dir.path,
                         column_specification=None)

            execution_context = ["##jacquard.version={0}".format(utils.__version__),
                "##jacquard.command=foo",
                "##jacquard.cwd=bar"]
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
                    self.assertTrue(actual[i] == "##jacquard.cwd=foo")
                elif expected[i].startswith("##jacquard.command="):
                    self.assertTrue(actual[i] == "##jacquard.command=bar")
                else:
                    self.assertEquals(expected[i], actual[i])

class ExpandFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_expand(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "06_expand")
            input_dir = os.path.join(module_testdir, "input")

            command = ["expand", input_dir, output_dir.path, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

    def test_expand_colSpec(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "06_expand_col_spec")
            input_dir = os.path.join(module_testdir, "input")
            col_spec = os.path.join(test_dir, "functional_tests", "col_spec.txt")
            command = ["expand", input_dir, output_dir.path, "--force", "--column_specification=" + col_spec]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)


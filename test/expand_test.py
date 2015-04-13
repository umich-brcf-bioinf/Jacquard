#pylint: disable=line-too-long, too-many-instance-attributes, unused-argument
#pylint: disable=global-statement, star-args, too-few-public-methods, invalid-name
#pylint: disable=too-many-public-methods
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import os

from testfixtures import TempDirectory

import jacquard.expand as expand
import jacquard.utils.logger
import jacquard.utils.utils as utils
import jacquard.utils.vcf as vcf
import test.utils.mock_logger
import test.utils.test_case as test_case
from test.utils.vcf_test import MockVcfReader


class ExpandTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(ExpandTestCase, self).setUp()
        expand.logger = test.utils.mock_logger

    def tearDown(self):
        test.utils.mock_logger.reset()
        expand.mock_logger = jacquard.utils.logger
        super(ExpandTestCase, self).tearDown()

    def test_create_row_dict(self):
        column_list = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                       "INFO", "FORMAT", "SAMPLE_A|NORMAL", "SAMPLE_A|TUMOR"]
        sample_tag_values = {"SAMPLE_A|NORMAL":{"DP":"50", "AF":"0.2"},
                             "SAMPLE_A|TUMOR":{"DP":"87", "AF":"0.3"}}
        vcf_record = vcf.VcfRecord("1", "42", "A", "AT",
                                   vcf_id="rs32", qual="30", vcf_filter="PASS",
                                   info="SNP;SOMATIC=1",
                                   sample_tag_values=sample_tag_values)
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

    def test_create_actual_column_list_regexImplictlyAnchored(self):
        potential_col_list = ["FOO", "FOO_1", "FOO_2", "BAR", "BAR_FOO"]
        col_spec = ["FOO"]
        actual_column_list = expand._create_actual_column_list(col_spec,
                                                               potential_col_list,
                                                               "col_spec.txt")
        expected_column_list = ["FOO"]

        self.assertEquals(expected_column_list, actual_column_list)


    def test_create_actual_column_list_duplicateRegexMatchDoesNotDuplicateColumns(self):
        potential_col_list = ["CHROM",
                              "JQ_DP|SAMPLE_A|NORMAL",
                              "JQ_DP|SAMPLE_B|NORMAL",
                              "JQ_AF|SAMPLE_A|NORMAL",
                              "JQ_AF|SAMPLE_B|NORMAL"]
        col_spec = ["CHROM", r"JQ_DP.*", "JQ_.*"]
        actual_column_list = expand._create_actual_column_list(col_spec,
                                                               potential_col_list,
                                                               "col_spec.txt")
        expected_column_list = ["CHROM",
                                "JQ_DP|SAMPLE_A|NORMAL",
                                "JQ_DP|SAMPLE_B|NORMAL",
                                "JQ_AF|SAMPLE_A|NORMAL",
                                "JQ_AF|SAMPLE_B|NORMAL"]

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

        actual_log_warnings = test.utils.mock_logger.messages["WARNING"]
        self.assertEquals([exp_warning_1, exp_warning_2],
                          actual_log_warnings)

    def test_create_potential_column_list(self):
        metaheaders = ['##INFO=<ID=AF,Number=1>',
                       '##INFO=<ID=AA,Number=1>',
                       '##FORMAT=<ID=GT,Number=1>',
                       '##FORMAT=<ID=GQ,Number=1,Description="bar">']
        sample_names = ["sampleA", "sampleB"]
        column_header = self.entab("#chrom|pos|id|ref|alt|"
                                   "qual|filter|info|format|"
                                   "sampleA|sampleB")
        mock_vcf_reader = MockVcfReader(metaheaders=metaheaders,
                                        column_header=column_header,
                                        sample_names=sample_names)

        actual_col_list = expand._create_potential_column_list(mock_vcf_reader)
        expected_col_list = ["chrom", "pos", "id", "ref", "alt",
                             "qual", "filter", "info",
                             "AA", "AF",
                             "GQ|sampleA", "GQ|sampleB",
                             "GT|sampleA", "GT|sampleB"]
        self.assertEquals(expected_col_list, actual_col_list)

    def test_create_potential_column_list_preservesSampleOrdering(self):
        metaheaders = ['##FORMAT=<ID=B,Number=1>', '##FORMAT=<ID=A,Number=1>']
        sample_names = ["sample1", "sample2", "sample10"]
        column_header = self.entab("#chrom|pos|id|ref|alt|"
                                   "qual|filter|info|format|"
                                   "sample1|sample2|sample10")

        mock_vcf_reader = MockVcfReader(metaheaders=metaheaders,
                                        column_header=column_header,
                                        sample_names=sample_names)

        actual_col_list = expand._create_potential_column_list(mock_vcf_reader)
        actual_format_sample_names = actual_col_list[8:]

        expected_format_sample_names = ["A|sample1", "A|sample2", "A|sample10",
                                        "B|sample1", "B|sample2", "B|sample10"]
        self.assertEquals(expected_format_sample_names,
                          actual_format_sample_names)

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

    def test_execute(self):
        vcf_content = ('''##source=strelka
##INFO=<ID=SOMATIC,Number=1,Description="foo">
##FORMAT=<ID=GT,Number=1,Description="bar">
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR
chr1|1|.|A|C|.|.|SOMATIC|GT|0/1|0/1
chr2|1|.|A|C|.|.|SOMATIC|GT|0/1|0/1
''').replace('|', "\t")

        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("P1.vcf", vcf_content.encode("utf8"))
            input_file = os.path.join(input_dir.path, "P1.vcf")
            output_file = os.path.join(output_dir.path, "P1.txt")
            args = Namespace(input=input_file,
                             original_output=output_file,
                             output=output_file,
                             column_specification=0)

            expand.execute(args, ["##extra_header1", "##extra_header2"])

            output_dir.check("P1.txt")
            actual_filename = os.path.join(output_dir.path, "P1.txt")
            with open(actual_filename) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()

        self.assertEquals(3, len(actual_output_lines))


    def test_predict_output(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("foo.txt", b"")
            args = Namespace(input=os.path.join(input_dir.path, "foo.txt"),
                             output=os.path.join(output_dir.path,
                                                 "expanded.txt"))

            desired_output_files = expand._predict_output(args)
            expected_desired_output_files = set(["expanded.txt"])

            self.assertEquals(expected_desired_output_files,
                              desired_output_files)

    def test_validate_args_colSpecValid(self):
        with TempDirectory() as col_spec_dir:
            col_spec_dir.write("col_spec.txt", b"chrom\npos\ninfo")
            col_spec_file = os.path.join(col_spec_dir.path, "col_spec.txt")

            args = Namespace(input="input.txt",
                             output="expanded.txt",
                             column_specification=col_spec_file)
            expand.validate_args(args)
            self.ok()

    def test_validate_args_colSpecInvalid(self):
        with TempDirectory() as col_spec_dir:
            col_spec_dir.write("col_spec.txt", b"chrom\npos\ninfo")

            args = Namespace(input="input.txt",
                             output="expanded.txt",
                             column_specification=col_spec_dir.path)
            self.assertRaisesRegexp(utils.UsageError,
                                    "The column specification file .* could not be read. Review inputs/usage and try again",
                                    expand.validate_args,
                                    args)

class ExpandFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_expand(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir,
                                          "functional_tests",
                                          "05_expand")
            input_dir = os.path.join(module_testdir, "input")

            command = ["expand",
                       os.path.join(input_dir, "summarized.vcf"),
                       os.path.join(output_dir.path, "expanded.tsv"),
                       "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")
            self.assertCommand(command, expected_dir)

    def test_expand_colSpec(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir,
                                          "functional_tests",
                                          "05_expand_col_spec")
            input_dir = os.path.join(module_testdir, "input")
            col_spec = os.path.join(test_dir, "functional_tests", "col_spec.txt")
            command = ["expand",
                       os.path.join(input_dir, "summarized.vcf"),
                       os.path.join(output_dir.path, "expanded.txt"),
                       "--force",
                       "--column_specification=" + col_spec]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

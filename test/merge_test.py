#pylint: disable=too-many-lines,missing-docstring,line-too-long,too-many-public-methods
#pylint: disable=too-few-public-methods,too-many-instance-attributes
#pylint: disable=too-many-arguments,invalid-name,protected-access,
#pylint: disable=global-statement
from __future__ import print_function, absolute_import, division

import argparse
from argparse import Namespace
from collections import OrderedDict
import os

from testfixtures import TempDirectory

import jacquard.utils.logger
import jacquard.merge as merge
import jacquard.utils.utils as utils
from jacquard.utils.vcf import VcfRecord
import jacquard.utils.vcf as vcf
import test.utils.mock_logger
import test.utils.test_case as test_case
from test.utils.vcf_test import MockVcfReader, MockFileReader, MockFileWriter, MockVcfRecord


class MockBufferedReader(object):
    def __init__(self, vcf_records):
        self.vcf_records_iter = iter(vcf_records)

    def next_if_equals(self, dummy):
        return next(self.vcf_records_iter)

class MergeVcfReaderTestCase(test_case.JacquardBaseTestCase):
    def test_extends_vcf_readers(self):
        file_contents = ["##metaheader1\n",
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|DP|32|78")]
        file_reader = MockFileReader("A.mutect.vcf", file_contents)
        vcf_reader = vcf.VcfReader(file_reader)
        merge_vcf_reader = merge.MergeVcfReader(vcf_reader._file_reader)

        self.assertEquals(merge_vcf_reader._file_reader, file_reader)

    def test_modify_metaheader(self):
        file_contents = ["##metaheader1\n",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        merge_vcf_reader = merge.MergeVcfReader(mock_file_reader)

        original_metaheader = '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">'
        transformed_tag = "JX1_DP"

        self.assertIn('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">', merge_vcf_reader.metaheaders)

        merge_vcf_reader.modify_metaheader(original_metaheader, transformed_tag)

        self.assertIn('##FORMAT=<ID=JX1_DP,Number=1,Type=Integer,Description="Read Depth">', merge_vcf_reader.metaheaders)
        self.assertNotIn('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">', merge_vcf_reader.metaheaders)

    def test_modify_format_tags(self):
        file_contents = ["##metaheader1\n",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        merge_vcf_reader = merge.MergeVcfReader(mock_file_reader)

        format_tags = {"DP": "JX1_DP"}
        vcf_record = MockVcfRecord("chr1", "245", "A", "G", vcf_format="AF:DP", samples=["0.2:21", "0.34:56"])

        self.assertEquals(OrderedDict({0: {"AF": "0.2", "DP": "21"}, 1: {"AF": "0.34", "DP": "56"}}),
                          vcf_record.sample_tag_values)

        merge_vcf_reader.modify_format_tag(vcf_record, format_tags)

        self.assertEquals(OrderedDict({0: {"AF": "0.2", "JX1_DP": "21"}, 1: {"AF": "0.34", "JX1_DP": "56"}}),
                          vcf_record.sample_tag_values)

    def test_modify_format_tags_inconsistentFormatTags(self):
        file_contents = ["##metaheader1\n",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|AF:DP|NORMAL|TUMOR"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|AF:DP:GT|NORMAL|TUMOR")]
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        merge_vcf_reader = merge.MergeVcfReader(mock_file_reader)

        format_tags = {"DP": "JX1_DP", "GT": "GT"}
        vcf_record = MockVcfRecord("chr1", "245", "A", "G", vcf_format="AF:DP", samples=["0.2:21", "0.34:56"])

        self.assertEquals(OrderedDict({0: {"AF": "0.2", "DP": "21"}, 1: {"AF": "0.34", "DP": "56"}}),
                          vcf_record.sample_tag_values)

        merge_vcf_reader.modify_format_tag(vcf_record, format_tags)

        self.assertEquals(OrderedDict({0: {"AF": "0.2", "JX1_DP": "21"}, 1: {"AF": "0.34", "JX1_DP": "56"}}),
                          vcf_record.sample_tag_values)

    def test_vcf_records_modifiesFormatTags(self):
        file_contents = ["##metaheader1\n",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|AF:DP|0.2:21|0.34:56")]
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        merge_vcf_reader = merge.MergeVcfReader(mock_file_reader)
        merge_vcf_reader.open()

        format_tags = {"DP": "JX1_DP"}

        vcf_records = merge_vcf_reader.vcf_records(format_tags, qualified=False)
        for vcf_record in vcf_records:
            normal_dict = OrderedDict(sorted({"AF": "0.2", "JX1_DP": "21"}.items()))
            tumor_dict = OrderedDict(sorted({"AF": "0.34", "JX1_DP": "56"}.items()))
            self.assertEquals(OrderedDict(sorted({"SampleNormal": normal_dict,
                                                  "SampleTumor": tumor_dict}.items())),
                              vcf_record.sample_tag_values)
        merge_vcf_reader.close()

    def test_store_unambiguous_format_tags(self):
        file_contents = ["##metaheader1\n",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|AF:DP|0.2:21|0.34:56")]
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        merge_vcf_reader = merge.MergeVcfReader(mock_file_reader)

        original_tag = "DP"
        new_tag = "JX1_DP"

        merge_vcf_reader.store_format_tags(original_tag, new_tag)
        self.assertEquals({"DP": "JX1_DP"}, merge_vcf_reader.format_tags)

class FilterTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(FilterTestCase, self).setUp()
        merge.logger = test.utils.mock_logger

    def tearDown(self):
        test.utils.mock_logger.reset()
        merge.logger = jacquard.utils.logger
        super(FilterTestCase, self).tearDown()

    def test_init_includeAllFlag(self):
        args = Namespace(include_all=True, include_cells=False, include_rows=False)
        record_filter = merge._Filter(args)
        self.assertEquals(merge._Filter._include_cell_if_all,
                          record_filter._cell_filter_strategy)
        self.assertEquals(merge._Filter._include_row_if_all,
                          record_filter._row_filter_strategy)

#     def test_init_includeAllFlag_raisesError(self):
#         args = Namespace(include_all=True, include_cells=True, include_rows=True)
#         self.assertRaisesRegexp(utils.UsageError,
#                                 "Unable to process command-line arguments. Neither --include_cells nor --include_rows can be specified if --include_all is specified.",
#                                  merge._Filter,
#                                  args)
#
#         args = Namespace(include_all=True, include_cells=False, include_rows=True)
#         self.assertRaisesRegexp(utils.UsageError,
#                                 "Unable to process command-line arguments. Neither --include_cells nor --include_rows can be specified if --include_all is specified.",
#                                  merge._Filter,
#                                  args)
#
#         args = Namespace(include_all=True, include_cells=True, include_rows=False)
#         self.assertRaisesRegexp(utils.UsageError,
#                                 "Unable to process command-line arguments. Neither --include_cells nor --include_rows can be specified if --include_all is specified.",
#                                  merge._Filter,
#                                  args)

    def test_init_includeValidAnysomaticByDefault(self):
        args = Namespace(include_all=False, include_cells="valid", include_rows="at_least_one_somatic")
        record_filter = merge._Filter(args)
        self.assertEquals(merge._Filter._include_cell_if_valid,
                          record_filter._cell_filter_strategy)
        self.assertEquals(merge._Filter._include_row_if_any_somatic,
                          record_filter._row_filter_strategy)

    def test_init_includeVariantPassed(self):
        args = Namespace(include_all=False, include_cells="passed", include_rows="all")
        record_filter = merge._Filter(args)
        self.assertEquals(merge._Filter._include_cell_if_passed,
                          record_filter._cell_filter_strategy)

    def test_init_includeVariantSomatic(self):
        args = Namespace(include_all=False, include_cells="somatic", include_rows="all")
        record_filter = merge._Filter(args)
        self.assertEquals(merge._Filter._include_cell_if_somatic,
                          record_filter._cell_filter_strategy)

    def test_init_includeLocusAllPassed(self):
        args = Namespace(include_all=False, include_cells="all", include_rows="all_passed")
        record_filter = merge._Filter(args)
        self.assertEquals(merge._Filter._include_row_if_all_passed,
                          record_filter._row_filter_strategy)

    def test_init_includeLocusAnyPassed(self):
        args = Namespace(include_all=False, include_cells="all", include_rows="at_least_one_passed")
        record_filter = merge._Filter(args)
        self.assertEquals(merge._Filter._include_row_if_any_passed,
                          record_filter._row_filter_strategy)

    def test_init_includeLocusAllSomatic(self):
        args = Namespace(include_all=False, include_cells="all", include_rows="all_somatic")
        record_filter = merge._Filter(args)
        self.assertEquals(merge._Filter._include_row_if_all_somatic,
                          record_filter._row_filter_strategy)

    def test_init_includeLocusAnySomatic(self):
        args = Namespace(include_all=False, include_cells="all", include_rows="at_least_one_somatic")
        record_filter = merge._Filter(args)
        self.assertEquals(merge._Filter._include_row_if_any_somatic,
                          record_filter._row_filter_strategy)

    def test_include_all(self):
        args = Namespace(include_all=False, include_rows = "all", include_cells = "all")
        filter_obj = merge._Filter(args)
        self.assertEquals(True, filter_obj._include_cell_if_all(None))
        self.assertEquals(True, filter_obj._include_row_if_all(None))

    def test_include_cell_if_valid(self):
        args = Namespace(include_all=False, include_rows = "all", include_cells = "all")
        filter_obj= merge._Filter(args)
        rec = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="JQ_EXCLUDE_FOO")
        self.assertEquals(False, filter_obj._include_cell_if_valid(rec))

        rec = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="bar")
        self.assertEquals(True, filter_obj._include_cell_if_valid(rec))

        self.assertEquals({"JQ_EXCLUDE_FOO": 1}, filter_obj.excluded_breakdown)

    def test_include_cell_if_passed(self):
        args = Namespace(include_all=False, include_cells = "all", include_rows = "all")
        filter_obj = merge._Filter(args)
        rec = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="foo")
        self.assertEquals(False, filter_obj._include_cell_if_passed(rec))

        rec = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="PASS")
        self.assertEquals(True, filter_obj._include_cell_if_passed(rec))

        self.assertEquals({"foo": 1}, filter_obj.excluded_breakdown)

    def test_include_cell_if_somatic(self):
        args = Namespace(include_all=False, include_rows = "all", include_cells = "all")
        filter_obj = merge._Filter(args)
        rec = VcfRecord("chrom", "pos", "ref", "alt",
                        sample_tag_values={"SA": {merge._JQ_SOMATIC_TAG: "0"}})
        self.assertEquals(False, filter_obj._include_cell_if_somatic(rec))

        rec = VcfRecord("chrom", "pos", "ref", "alt",
                        sample_tag_values={"SA": {merge._JQ_SOMATIC_TAG: "1"}})
        self.assertEquals(True, filter_obj._include_cell_if_somatic(rec))

        self.assertEquals({"not somatic": 1}, filter_obj.excluded_breakdown)

    def test_include_row_if_all_passed(self):
        rec1 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="PASS")
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="PASS")
        self.assertEquals(True,
                          merge._Filter._include_row_if_all_passed([rec1, rec2]))

        rec1 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="PASS")
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="FAIL")
        self.assertEquals(False,
                          merge._Filter._include_row_if_all_passed([rec1, rec2]))

        rec1 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="FAIL")
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="FAIL")
        self.assertEquals(False,
                          merge._Filter._include_row_if_all_passed([rec1, rec2]))

    def test_include_row_if_any_passed(self):
        rec1 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="PASS")
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="PASS")
        self.assertEquals(True,
                          merge._Filter._include_row_if_any_passed([rec1, rec2]))

        rec1 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="PASS")
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="FAIL")
        self.assertEquals(True,
                          merge._Filter._include_row_if_any_passed([rec1, rec2]))

        rec1 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="FAIL")
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="FAIL")
        self.assertEquals(False,
                          merge._Filter._include_row_if_any_passed([rec1, rec2]))

    def test_include_row_if_all_somatic(self):
        nonsomatic = {"SA": {merge._JQ_SOMATIC_TAG: "0"}}
        somatic = {"SA": {merge._JQ_SOMATIC_TAG: "1"}}
        rec1 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=somatic)
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=somatic)
        self.assertEquals(True,
                          merge._Filter._include_row_if_all_somatic([rec1, rec2]))

        rec1 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=somatic)
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=nonsomatic)
        self.assertEquals(False,
                          merge._Filter._include_row_if_all_somatic([rec1, rec2]))

        rec1 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=nonsomatic)
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=nonsomatic)
        self.assertEquals(False,
                          merge._Filter._include_row_if_all_somatic([rec1, rec2]))

    def test_include_row_if_any_somatic(self):
        nonsomatic = {"SA": {merge._JQ_SOMATIC_TAG: "0"}}
        somatic = {"SA": {merge._JQ_SOMATIC_TAG: "1"}}
        rec1 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=somatic)
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=somatic)
        self.assertEquals(True,
                          merge._Filter._include_row_if_any_somatic([rec1, rec2]))

        rec1 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=somatic)
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=nonsomatic)
        self.assertEquals(True,
                          merge._Filter._include_row_if_any_somatic([rec1, rec2]))

        rec1 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=nonsomatic)
        rec2 = VcfRecord("chrom", "pos", "ref", "alt", sample_tag_values=nonsomatic)
        self.assertEquals(False,
                          merge._Filter._include_row_if_any_somatic([rec1, rec2]))

    def test_include_row_statistics(self):
        args = Namespace(include_all=False,include_cells="all", include_rows="at_least_one_passed")
        filter_obj = merge._Filter(args)
        passed_record = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="PASS")
        failed_record = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="FAIL")

        self.assertEquals(0, filter_obj.row_count)
        self.assertEquals(0, filter_obj.rows_excluded)

        filter_obj.include_row([passed_record, failed_record])
        self.assertEquals(1, filter_obj.row_count)
        self.assertEquals(0, filter_obj.rows_excluded)

        filter_obj.include_row([failed_record, failed_record])
        self.assertEquals(2, filter_obj.row_count)
        self.assertEquals(1, filter_obj.rows_excluded)

    def test_include_cell_statistics(self):
        args = Namespace(include_all=False,include_cells="passed", include_rows="all")
        filter_obj = merge._Filter(args)
        passed_record = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="PASS")
        failed_record = VcfRecord("chrom", "pos", "ref", "alt", vcf_filter="FAIL")

        self.assertEquals(0, filter_obj.cell_count)
        self.assertEquals(0, filter_obj.cells_excluded)

        filter_obj.include_cell(passed_record)
        self.assertEquals(1, filter_obj.cell_count)
        self.assertEquals(0, filter_obj.cells_excluded)

        filter_obj.include_cell(failed_record)
        self.assertEquals(2, filter_obj.cell_count)
        self.assertEquals(1, filter_obj.cells_excluded)

    def test_log_statistics(self):
        args = Namespace(include_all=False,include_cells="passed", include_rows="at_least_one_somatic")
        filter_obj = merge._Filter(args)
        filter_obj.cell_count = 20
        filter_obj.cells_excluded = 4
        filter_obj.row_count = 40
        filter_obj.rows_excluded = 2
        filter_obj.excluded_breakdown = {"JQ_EXCLUDED": 2, "JQ_INVALID": 1, "FOO": 1}

        filter_obj.log_statistics()
        actual_info_message = test.utils.mock_logger.messages["INFO"]
        actual_debug_message = test.utils.mock_logger.messages["DEBUG"]

        self.assertEquals(2, len(actual_info_message))
        self.assertEquals("20% (4) cells were excluded because (--include_cells=passed)",
                          actual_info_message[0])
        self.assertEquals("5% (2) rows were excluded because (--include_rows=at_least_one_somatic)",
                          actual_info_message[1])

        self.assertEquals(3, len(actual_debug_message))
        self.assertIn("2 cells were excluded with [JQ_EXCLUDED]",
                          actual_debug_message)
        self.assertIn("1 cells were excluded with [JQ_INVALID]",
                          actual_debug_message)
        self.assertIn("1 cells were excluded with [FOO]",
                          actual_debug_message)

class MergeTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(MergeTestCase, self).setUp()
        merge.logger = test.utils.mock_logger

    def tearDown(self):
        test.utils.mock_logger.reset()
        merge.logger = jacquard.utils.logger
        super(MergeTestCase, self).tearDown()

    def test_validate_consistent_arguments(self):
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(title="subcommands",
                                           dest="subparser_name")
        merge.add_subparser(subparsers)

        cell_choices = subparsers.choices["merge"]._option_string_actions["--include_cells"].choices
        row_choices = subparsers.choices["merge"]._option_string_actions["--include_rows"].choices

        cell_help = subparsers.choices["merge"]._option_string_actions["--include_cells"].help
        cell_help = [i.split(":")[0] for i in cell_help.split("\n")]
        row_help = subparsers.choices["merge"]._option_string_actions["--include_rows"].help
        row_help = [i.split(":")[0] for i in row_help.split("\n")]

        args = Namespace(include_all=False, include_rows="all", include_cells="all")
        filter_obj = merge._Filter(args)

        for key in filter_obj._cell_filters.keys():
            if key not in cell_choices:
                self.assertFalse("[{}] not found in {}".format(key, cell_choices))
            if key not in cell_help:
                self.assertFalse("[{}] not found in {}".format(key, cell_help))

        for key in filter_obj._row_filters.keys():
            if key not in row_choices:
                self.assertFalse("[{}] not found in {}".format(key, row_choices))
            if key not in row_help:
                self.assertFalse("[{}] not found in {}".format(key, row_help))

    def test_validate_consistent_samples_missingCaller(self):
        input_files = [MockFileReader("A.mutect.vcf",
                                      ["##jacquard.translate.caller=MuTect"]),
                       MockFileReader("A.strelka.vcf",
                                      ["##jacquard.translate.caller=Strelka"]),
                       MockFileReader("A.varscan.vcf",
                                      ["##jacquard.translate.caller=VarScan"]),
                       MockFileReader("B.strelka.vcf",
                                      ["##jacquard.translate.caller=Strelka"]),
                       MockFileReader("B.varscan.vcf",
                                      ["##jacquard.translate.caller=VarScan"]),
                       MockFileReader("C.varscan.vcf",
                                      ["##jacquard.translate.caller=VarScan"])]
        merge._validate_consistent_samples(input_files)

        actual_log_warnings = test.utils.mock_logger.messages["WARNING"]
        expected_log_warnings = ["Sample [B] is missing VCF(s): ['MuTect']",
                                 "Sample [C] is missing VCF(s): ['MuTect', 'Strelka']",
                                 "Some samples appear to be missing VCF(s)"]
        self.assertEquals(expected_log_warnings, actual_log_warnings)

    def test_validate_consistent_samples_allMissingCallers(self):
        input_files = [MockFileReader("A.mutect.vcf",
                                      ["##jacquard.translate.caller=MuTect"]),
                       MockFileReader("B.strelka.vcf",
                                      ["##jacquard.translate.caller=Strelka"]),
                       MockFileReader("C.mutect.vcf",
                                      ["##jacquard.translate.caller=MuTect"])]
        merge._validate_consistent_samples(input_files)

        actual_log_warnings = test.utils.mock_logger.messages["WARNING"]
        expected_log_warnings = ["Sample [A] is missing VCF(s): ['Strelka']",
                                 "Sample [B] is missing VCF(s): ['MuTect']",
                                 "Sample [C] is missing VCF(s): ['Strelka']",
                                 "Some samples appear to be missing VCF(s)"]
        self.assertEquals(expected_log_warnings, actual_log_warnings)

    def test_validate_consistent_input_allTranslatedOkay(self):
        metaheaders1 = ["##metaheader1", "##metaheader2", "##jacquard.translate.caller"]
        vcf_reader1 = MockVcfReader(metaheaders=metaheaders1)

        metaheaders2 = ["##metaheader3", "##jacquard.translate.caller", "##metaheader4"]
        vcf_reader2 = MockVcfReader(metaheaders=metaheaders2)

        merge._validate_consistent_input([vcf_reader1, vcf_reader2], False)
        actual_log_errors = test.utils.mock_logger.messages["ERROR"]
        self.assertEquals([], actual_log_errors)

    def test_validate_consistent_input_allUntranslatedOkay(self):
        metaheaders1 = ["##metaheader1", "##metaheader2"]
        vcf_reader1 = MockVcfReader(metaheaders=metaheaders1)

        metaheaders2 = ["##metaheader3", "##metaheader4"]
        vcf_reader2 = MockVcfReader(metaheaders=metaheaders2)

        merge._validate_consistent_input([vcf_reader1, vcf_reader2], False)
        actual_log_errors = test.utils.mock_logger.messages["ERROR"]
        self.assertEquals([], actual_log_errors)

    def test_validate_consistent_input_errorIfMixed(self):
        metaheaders1 = ["##metaheader1", "##metaheader2", "##jacquard.translate.caller"]
        vcf_reader1 = MockVcfReader(metaheaders=metaheaders1)

        metaheaders2 = ["##metaheader3", "##metaheader4"]
        vcf_reader2 = MockVcfReader(metaheaders=metaheaders2)

        self.assertRaisesRegexp(utils.UsageError,
                                r"Some input VCFs \[.*\] were not translated by Jacquard. Review input and/or use --included_all flag",
                                merge._validate_consistent_input,
                                [vcf_reader1, vcf_reader2],
                                False)

    def test_validate_consistent_input_okayIfIncludeAll(self):
        metaheaders1 = ["##metaheader1", "##metaheader2", "##jacquard.translate.caller"]
        vcf_reader1 = MockVcfReader(metaheaders=metaheaders1)

        metaheaders2 = ["##metaheader3", "##metaheader4"]
        vcf_reader2 = MockVcfReader(metaheaders=metaheaders2)

        merge._validate_consistent_input([vcf_reader1, vcf_reader2], True)
        actual_log_errors = test.utils.mock_logger.messages["ERROR"]
        self.assertEquals([], actual_log_errors)


    def test_predict_output(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.normalized.jacquardTags.HCsomatic.vcf", b"##source=strelka\n#colHeader")
            input_dir.write("B.normalized.jacquardTags.HCsomatic.vcf", b"##source=strelka\n#colHeader")
            output_dir.write("merged.vcf", b"##source=strelka\n#colHeader")
            args = Namespace(input=input_dir.path,
                             output=os.path.join(output_dir.path, "merged.vcf"))

            desired_output_files = merge._predict_output(args)
            expected_desired_output_files = set(["merged.vcf"])

            self.assertEquals(expected_desired_output_files, desired_output_files)

    def test_get_caller_name(self):
        file_contents = ["##metaheader1\n",
                         "##jacquard.translate.caller=MuTect",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        merge_vcf_reader = merge.MergeVcfReader(mock_file_reader)

        caller_name = merge._get_caller_name(merge_vcf_reader)

        self.assertEquals("MuTect", caller_name)

    def test_get_caller_name_noCallerDefined(self):
        file_contents = ["##metaheader1\n",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        merge_vcf_reader = merge.MergeVcfReader(mock_file_reader)

        caller_name = merge._get_caller_name(merge_vcf_reader)

        self.assertEquals(False, caller_name)

    def test_alter_description(self):
        metaheadr = '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'
        caller_name = 'VarScan'

        expected = '##FORMAT=<ID=AF,Number=A,Type=Float,Description="[VarScan]: Allele Frequency">'
        actual = merge._alter_description(metaheadr, caller_name)
        self.assertEquals(expected, actual)

    def test_get_format_tag_regex_defaultValue(self):
        args = Namespace(include_all=False, tags=False)
        actual_regex = merge._get_format_tag_regex(args)

        self.assertEquals(merge._DEFAULT_INCLUDED_FORMAT_TAGS, actual_regex)

    def test_get_format_tag_regex_specifiedValue(self):
        args = Namespace(include_all=False, tags="^FOO.*,BAR.*")
        actual_regex = merge._get_format_tag_regex(args)

        self.assertEquals(["^FOO.*", "BAR.*"], actual_regex)

    def test_get_format_tag_regex_includeAll(self):
        args = Namespace(include_all=True, tags=False)
        actual_regex = merge._get_format_tag_regex(args)

        self.assertEquals([".*"], actual_regex)

    def test_get_format_tag_regex_raisesError(self):
        args = Namespace(include_all=True, tags=True)
        self.assertRaisesRegexp(utils.UsageError,
                                "Unable to process command-line arguments. --include_format_tags cannot be specified if --include_all is specified.",
                                merge._get_format_tag_regex,
                                args)

    def test_get_format_tags(self):
        file_contents = ["##metaheader1\n",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_file_reader = MockFileReader("my_dir/my_file.txt", file_contents)
        merge_vcf_reader = merge.MergeVcfReader(mock_file_reader)

        format_tags = merge._get_format_tags([merge_vcf_reader])
        expected_format_tag_dict = {"AF": ['##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'],
                                    "DP": ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">']}

        self.assertEquals(expected_format_tag_dict, format_tags)

    def test_get_format_tags_considersDifferentCallers(self):
        file_contents1 = ["##metaheader1\n",
                          "##jacquard.translate.caller=MuTect",
                          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                          self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                          self.entab("chr2|1|.|A|C|.|.|INFO|DP|15|87")]
        mock_file_reader1 = MockFileReader("fileA.txt", file_contents1)
        merge_vcf_reader1 = merge.MergeVcfReader(mock_file_reader1)

        file_contents2 = ["##metaheader1\n",
                          "##jacquard.translate.caller=VarScan",
                          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                          self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                          self.entab("chr2|1|.|A|C|.|.|INFO|DP|75|65")]
        mock_file_reader2 = MockFileReader("fileB.txt", file_contents2)
        merge_vcf_reader2 = merge.MergeVcfReader(mock_file_reader2)

        format_tags = merge._get_format_tags([merge_vcf_reader1, merge_vcf_reader2])
        expected_format_tag_dict = {"DP": ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="[MuTect]: Read Depth">',
                                           '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="[VarScan]: Read Depth">']}

        self.assertEquals(expected_format_tag_dict, format_tags)

    def test_get_format_tags_sameMetaheaderOkay(self):
        file_contents1 = ["##metaheader1\n",
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_file_reader1 = MockFileReader("my_dir/my_file.txt", file_contents1)
        merge_vcf_reader1 = merge.MergeVcfReader(mock_file_reader1)

        file_contents2 = ["##metaheader1\n",
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_file_reader2 = MockFileReader("my_dir/my_file.txt", file_contents2)
        merge_vcf_reader2 = merge.MergeVcfReader(mock_file_reader2)

        format_tags = merge._get_format_tags([merge_vcf_reader1, merge_vcf_reader2])
        expected_format_tag_dict = {"DP": ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">']}

        self.assertEquals(expected_format_tag_dict, format_tags)

    def test_disambiguate_format_tags_modifiesMetaheader(self):
        file_contents1 = ["##metaheader1\n",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_file_reader1 = MockFileReader("fileA.txt", file_contents1)

        merge_vcf_reader1 = merge.MergeVcfReader(mock_file_reader1)
        file_contents2 = ["##metaheader1\n",
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|FORMAT|NORMAL|TUMOR")]
        mock_file_reader2 = MockFileReader("fileB.txt", file_contents2)
        merge_vcf_reader2 = merge.MergeVcfReader(mock_file_reader2)

        format_tags = {"DP": ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate Read Depth">'],
                       "AF": ['##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">']}
        merge._disambiguate_format_tags([merge_vcf_reader1, merge_vcf_reader2], format_tags)

        self.assertIn('##FORMAT=<ID=JX1_DP,Number=1,Type=Integer,Description="Read Depth">',
                      merge_vcf_reader1.metaheaders)
        self.assertIn('##FORMAT=<ID=JX2_DP,Number=1,Type=Integer,Description="Approximate Read Depth">',
                      merge_vcf_reader2.metaheaders)

    def xtest_disambiguate_format_tags_considersDifferentCallers(self):
        file_contents1 = ["##metaheader1\n",
                          "##jacquard.translate.caller=MuTect",
                          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                          self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                          self.entab("chr2|1|.|A|C|.|.|INFO|DP|15|87")]
        mock_file_reader1 = MockFileReader("fileA.txt", file_contents1)
        merge_vcf_reader1 = merge.MergeVcfReader(mock_file_reader1)

        file_contents2 = ["##metaheader1\n",
                          "##jacquard.translate.caller=VarScan",
                          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                          self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                          self.entab("chr2|1|.|A|C|.|.|INFO|DP|75|65")]
        mock_file_reader2 = MockFileReader("fileB.txt", file_contents2)
        merge_vcf_reader2 = merge.MergeVcfReader(mock_file_reader2)

        format_tags = {"DP": ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate Read Depth">']}
        merge._disambiguate_format_tags([merge_vcf_reader1, merge_vcf_reader2], format_tags)

        self.assertIn('##FORMAT=<ID=JX1_DP,Number=1,Type=Integer,Description="[MuTect]: Read Depth">',
                      merge_vcf_reader1.metaheaders)
        self.assertIn('##FORMAT=<ID=JX2_DP,Number=1,Type=Integer,Description="[VarScan]: Read Depth">',
                      merge_vcf_reader2.metaheaders)

    def test_disambiguate_format_tags_modifiesRecords(self):
        file_contents1 = ["##metaheader1\n",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|Normal|Tumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|AF:DP|0.2:34|0.1:54")]
        mock_file_reader1 = MockFileReader("fileA.txt", file_contents1)

        merge_vcf_reader1 = merge.MergeVcfReader(mock_file_reader1)
        file_contents2 = ["##metaheader1\n",
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|Normal|Tumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|DP|12|76")]
        mock_file_reader2 = MockFileReader("fileB.txt", file_contents2)
        merge_vcf_reader2 = merge.MergeVcfReader(mock_file_reader2)

        format_tags = {"DP": ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate Read Depth">'],
                       "AF": ['##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">']}
        merge._disambiguate_format_tags([merge_vcf_reader1, merge_vcf_reader2], format_tags)

        self.assertEquals({"DP": "JX1_DP", "AF": "AF"}, merge_vcf_reader1.format_tags)
        self.assertEquals({"DP": "JX2_DP"}, merge_vcf_reader2.format_tags)


    def test_build_coordinates(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "G", "id=2")
        fileBrec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileArec1, fileArec2]),
                        MockVcfReader(records=[fileBrec1, fileBrec2])]

        actual_coordinates = merge._build_coordinates(mock_readers)

        expected = [fileArec1, fileArec2, fileBrec2]
        self.assertEquals(expected, actual_coordinates)

    def test_build_coordinates_multAltsEmpty(self):
        fileArec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileArec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileBrec1 = vcf.VcfRecord("chr2", "12", "A", "G", "id=2")
        fileBrec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileArec1, fileArec2]),
                        MockVcfReader(records=[fileBrec1, fileBrec2])]

        actual_coordinates = merge._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        self.assertEquals([], actual_multalts)

    def test_build_coordinates_flagsMultAltsFromDistinctFiles(self):
        fileA_rec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileA_rec2 = vcf.VcfRecord("chr2", "12", "A", "G", "id=1")
        fileB_rec1 = vcf.VcfRecord("chr2", "12", "A", "T", "id=2")
        fileB_rec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileA_rec1, fileA_rec2]),
                        MockVcfReader(records=[fileB_rec1, fileB_rec2])]

        actual_coordinates = merge._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = [fileA_rec2, fileB_rec1]
        self.assertEquals(expected, actual_multalts)

    def test_build_coordinates_flagsMultAltsWithinFile(self):
        fileA_rec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileA_rec2 = vcf.VcfRecord("chr2", "12", "A", "G,T", "id=1")
        fileB_rec1 = vcf.VcfRecord("chr3", "12", "A", "T", "id=2")
        fileB_rec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileA_rec1, fileA_rec2]),
                        MockVcfReader(records=[fileB_rec1, fileB_rec2])]

        actual_coordinates = merge._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = [fileA_rec2]
        self.assertEquals(expected, actual_multalts)

    def test_build_coordinates_flagsMultAltsWithDistinctRefs(self):
        fileA_rec1 = vcf.VcfRecord("chr1", "1", "A", "C")
        fileA_rec2 = vcf.VcfRecord("chr2", "2", "A", "G", "id=1")
        fileB_rec1 = vcf.VcfRecord("chr2", "2", "AT", "T", "id=2")
        fileB_rec2 = vcf.VcfRecord("chr42", "16", "G", "C")

        mock_readers = [MockVcfReader(records=[fileA_rec1, fileA_rec2]),
                        MockVcfReader(records=[fileB_rec1, fileB_rec2])]

        actual_coordinates = merge._build_coordinates(mock_readers)

        actual_multalts = [record for record in actual_coordinates if record.info == "JQ_MULT_ALT_LOCUS"]

        expected = [fileA_rec2, fileB_rec1]
        self.assertEquals(expected, actual_multalts)

    def test_build_merged_record_onlyKeepJQTags(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SD": {"JQ_foo":"bar1",
                              "blahJQ_": "bar2"},
                       "SC": {"JQ_foo":"bar3",
                              "blah":"bar4"}})
        samples2 = OD({"SB": {"JQ_foo":"bar5"},
                       "SA": {"JQ_foo":"bar6"}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        sample_list = ["SA", "SB", "SC", "SD"]
        tags_to_keep = ["JQ_foo", "JQ_foo"]
        actual_record = merge._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

        self.assertEquals(OD([("JQ_foo", "bar6")]), actual_record.sample_tag_values["SA"])
        self.assertEquals(OD([("JQ_foo", "bar5")]), actual_record.sample_tag_values["SB"])
        self.assertEquals(OD([("JQ_foo", "bar3")]), actual_record.sample_tag_values["SC"])
        self.assertEquals(OD([("JQ_foo", "bar1")]), actual_record.sample_tag_values["SD"])

    def test_build_merged_record_redundantPatientNames(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"PA|NORMAL": {"JQ_vs":"1"},
                       "PA|TUMOR": {"JQ_vs":"2"}})
        samples2 = OD({"PA|NORMAL": {"JQ_mt":"3"},
                       "PA|TUMOR": {"JQ_mt":"4"}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        sample_list = ["PA|NORMAL", "PA|TUMOR"]
        tags_to_keep = {"JQ_vs": "JQ_vs", "JQ_mt": "JQ_mt"}
        actual_record = merge._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

        self.assertEquals(OD([("JQ_mt", "3"), ("JQ_vs", "1")]), actual_record.sample_tag_values["PA|NORMAL"])
        self.assertEquals(OD([("JQ_mt", "4"), ("JQ_vs", "2")]), actual_record.sample_tag_values["PA|TUMOR"])

    def test_build_merged_record_preserveSampleNamesAndOrder(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SD": {"foo":"bar1"},
                       "SC": {"foo":"bar2"}})
        samples2 = OD({"SB": {"foo":"bar3"},
                       "SA": {"foo":"bar4"}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        sample_list = ["SD", "SA", "SC", "SB"]
        tags_to_keep = ["foo"]
        actual_record = merge._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

        self.assertEquals(sample_list, list(actual_record.sample_tag_values.keys()))

    def test_build_merged_record_fillsMissingSamples(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SA": {"foo":"bar3"},
                       "SB": {"foo":"bar4"}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)

        sample_list = ["SA", "SB", "SC", "SD"]
        tags_to_keep = ["foo"]
        actual_record = merge._build_merged_record(coordinate, [record1], sample_list, tags_to_keep)

        self.assertEquals(sample_list, list(actual_record.sample_tag_values.keys()))

    def test_build_merged_record_baseInfoCopiedFromCoordinate(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        samples1 = OD({"SA": {},
                       "SB": {}})
        samples2 = OD({"SC": {},
                       "SD": {}})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        actual_record = merge._build_merged_record(coordinate, [record1, record2], [], [])

        self.assertEquals("chr1", actual_record.chrom)
        self.assertEquals("1", actual_record.pos)
        self.assertEquals("A", actual_record.ref)
        self.assertEquals("C", actual_record.alt)
        self.assertEquals("baseInfo", actual_record.info)

    def test_build_merged_record_tagsOrdered(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        sampleA_tag_values = OD({"foo":"A1", "bar":"A2"})
        sampleB_tag_values = OD({"foo":"B1", "bar":"B2"})
        sampleC_tag_values = OD({"foo":"C1", "bar":"C2"})
        sampleD_tag_values = OD({"foo":"D1", "bar":"D2"})
        samples1 = OD({"SA": sampleA_tag_values,
                       "SB": sampleB_tag_values})
        samples2 = OD({"SC": sampleC_tag_values,
                       "SD": sampleD_tag_values})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        sample_list = ["SA", "SB", "SC", "SD"]
        tags_to_keep = ["foo", "bar"]
        actual_record = merge._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

        self.assertEquals(OD([("bar", "A2"), ("foo", "A1")]), actual_record.sample_tag_values["SA"])
        self.assertEquals(OD([("bar", "B2"), ("foo", "B1")]), actual_record.sample_tag_values["SB"])
        self.assertEquals(OD([("bar", "C2"), ("foo", "C1")]), actual_record.sample_tag_values["SC"])
        self.assertEquals(OD([("bar", "D2"), ("foo", "D1")]), actual_record.sample_tag_values["SD"])

    def test_build_merged_record_heterogeneousTags(self):
        OD = OrderedDict
        coordinate = VcfRecord("chr1", "1", "A", "C", info="baseInfo")
        sampleA_tag_values = OD({"foo":"A1", "bar":"A2"})
        sampleB_tag_values = OD({"foo":"B1", "bar":"B2"})
        sampleC_tag_values = OD({"baz":"C1"})
        sampleD_tag_values = OD({"baz":"D1"})
        samples1 = OD({"SA": sampleA_tag_values,
                       "SB": sampleB_tag_values})
        samples2 = OD({"SC": sampleC_tag_values,
                       "SD": sampleD_tag_values})
        record1 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples1)
        record2 = VcfRecord("chr1", "1", "A", "C", sample_tag_values=samples2)

        sample_list = ["SA", "SB", "SC", "SD"]
        tags_to_keep = ["foo", "bar", "baz"]
        actual_record = merge._build_merged_record(coordinate, [record1, record2], sample_list, tags_to_keep)

        self.assertEquals("chr1", actual_record.chrom)
        self.assertEquals("1", actual_record.pos)
        self.assertEquals("A", actual_record.ref)
        self.assertEquals("C", actual_record.alt)
        self.assertEquals("baseInfo", actual_record.info)
        self.assertEquals(set(["SA", "SB", "SC", "SD"]), set(actual_record.sample_tag_values.keys()))
        self.assertEquals(OD([("bar", "A2"), ("baz", "."), ("foo", "A1")]), actual_record.sample_tag_values["SA"])
        self.assertEquals(OD([("bar", "B2"), ("baz", "."), ("foo", "B1")]), actual_record.sample_tag_values["SB"])
        self.assertEquals(OD([("bar", "."), ("baz", "C1"), ("foo", ".")]), actual_record.sample_tag_values["SC"])
        self.assertEquals(OD([("bar", "."), ("baz", "D1"), ("foo", ".")]), actual_record.sample_tag_values["SD"])

    def test_merge_records(self):
        coordinates = [VcfRecord("chrom", "pos", "ref", "alt")]
        filter_strategy = merge._Filter(Namespace(include_all=False,
                                                  include_cells="all",
                                                  include_rows="all"))
        OD = OrderedDict
        record1 = VcfRecord("chrom", "pos", "ref", "alt",
                            sample_tag_values=OD({"SA": OD({"foo":"A"}),
                                                  "SB":OD({"foo":"B"})}))
        record2 = VcfRecord("chrom", "pos", "ref", "alt",
                            sample_tag_values=OD({"SC": OD({"foo":"C"}),
                                                  "SD":OD({"foo":"D"})}))
        vcf_readers = [MockVcfReader(input_filepath="fileA.vcf",
                                     records=[record1]),
                       MockVcfReader(input_filepath="fileB.vcf",
                                     records=[record2])]
        all_sample_names = ["SA", "SB", "SC", "SD"]
        format_tags_to_keep = ["foo"]
        file_writer = MockFileWriter()

        merge._merge_records(vcf_readers,
                             coordinates,
                             filter_strategy,
                             all_sample_names,
                             format_tags_to_keep,
                             file_writer)

        expected_record = "chrom\tpos\t.\tref\talt\t.\t.\t.\tfoo\tA\tB\tC\tD"
        self.assertEquals([expected_record], file_writer.lines())

    def test_pull_matching_records(self):
        filter_strategy = merge._Filter(Namespace(include_all=False,
                                                  include_cells="all",
                                                  include_rows="all"))
        coordinate = VcfRecord("chrom", "pos", "ref", "alt")
        OD = OrderedDict
        record1 = VcfRecord("chrom", "pos", "ref", "alt",
                            sample_tag_values=OD({"SA": OD({"foo":"A"}),
                                                  "SB":OD({"foo":"B"})}))
        record2 = VcfRecord("chrom", "pos", "ref", "alt",
                            sample_tag_values=OD({"SC": OD({"foo":"C"}),
                                                  "SD":OD({"foo":"D"})}))
        buffered_readers = [MockBufferedReader([record1]),
                            MockBufferedReader([record2])]

        vcf_records = merge._pull_matching_records(filter_strategy,
                                                   coordinate,
                                                   buffered_readers)
        self.assertEqual([record1, record2], vcf_records)

    def test_build_sample_list_simpleSampleList(self):
        reader1 = MockVcfReader("PA.foo.vcf",
                                sample_names=["Sample_A", "Sample_B"])
        reader2 = MockVcfReader("PA.bar.vcf",
                                sample_names=["Sample_A", "Sample_B"])
        reader3 = MockVcfReader("PB.vcf",
                                sample_names=["Sample_C", "Sample_D"])
        readers = [reader1, reader2, reader3]
        actual_sample_names, dummy = merge._build_sample_list(readers)

        expected_sample_names = ["PA|Sample_A",
                                 "PA|Sample_B",
                                 "PB|Sample_C",
                                 "PB|Sample_D"]
        self.assertEquals(expected_sample_names, actual_sample_names)

    def test_build_sample_list_patientNamesNaturalOrdered(self):
        reader1 = MockVcfReader("P10.foo.vcf", sample_names=["S10", "S2"])
        reader2 = MockVcfReader("P1A.bar.vcf", sample_names=["S10", "S2"])
        reader3 = MockVcfReader("P1.bar.vcf", sample_names=["S10", "S2"])
        reader4 = MockVcfReader("P2.vcf", sample_names=["S10", "S2"])
        readers = [reader1, reader2, reader3, reader4]
        actual_sample_names, dummy = merge._build_sample_list(readers)

        expected_sample_names = ["P1A|S2", "P1A|S10",
                                 "P1|S2", "P1|S10",
                                 "P2|S2", "P2|S10",
                                 "P10|S2", "P10|S10"]
        self.assertEquals(expected_sample_names, actual_sample_names)

    def test_build_sample_list_sampleNamesOrdered(self):
        reader1 = MockVcfReader("P10.foo.vcf", sample_names=["S10", "S1", "S100", "S2"])
        reader2 = MockVcfReader("P1.foo.vcf", sample_names=["S0", "S1"])
        readers = [reader1, reader2]
        actual_sample_names, dummy = merge._build_sample_list(readers)

        expected_sample_names = ["P1|S0", "P1|S1",
                                 "P10|S1", "P10|S2", "P10|S10", "P10|S100"]
        self.assertEquals(expected_sample_names, actual_sample_names)

    def test_build_sample_list_mergeMetaheaders(self):
        reader1 = MockVcfReader("PA.foo.vcf",
                                sample_names=["Sample_A", "Sample_B"])
        reader2 = MockVcfReader("PA.bar.vcf",
                                sample_names=["Sample_A", "Sample_B"])
        reader3 = MockVcfReader("PB.vcf",
                                sample_names=["Sample_C", "Sample_D"])
        readers = [reader1, reader2, reader3]

        dummy, actual_metaheaders = merge._build_sample_list(readers)

        expected_metaheaders = ["##jacquard.merge.sample=<Column=1,Name=PA|Sample_A,Source=PA.foo.vcf|PA.bar.vcf>",
                                "##jacquard.merge.sample=<Column=2,Name=PA|Sample_B,Source=PA.foo.vcf|PA.bar.vcf>",
                                "##jacquard.merge.sample=<Column=3,Name=PB|Sample_C,Source=PB.vcf>",
                                "##jacquard.merge.sample=<Column=4,Name=PB|Sample_D,Source=PB.vcf>"]
        self.assertEquals(4, len(actual_metaheaders))
        self.assertEquals(expected_metaheaders, actual_metaheaders)

    def test_create_vcf_readers(self):
        with TempDirectory() as input_dir:
            fileA = input_dir.write("fileA.vcf",
                                    (b"##source=strelka\n"
                                     b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_A\tSample_B\n"
                                     b"chr1\t31\t.\tA\tT\t.\t.\t.\tDP\t23\t52\n"))
            fileB = input_dir.write("fileB.vcf",
                                    (b"##source=strelka\n"
                                     b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_C\tSample_D\n"
                                     b"chr2\t32\t.\tA\tT\t.\t.\t.\tDP\t24\t53\n"))
            input_files = [vcf.FileReader(fileA), vcf.FileReader(fileB)]
            vcf_readers = merge._create_vcf_readers(input_files)

            for vcf_reader in vcf_readers:
                vcf_reader.close()

            self.assertEquals(2, len(vcf_readers))
            self.assertEquals("fileA.vcf", vcf_readers[0].file_name)
            self.assertEquals("fileB.vcf", vcf_readers[1].file_name)

    def test_create_buffered_readers(self):
        vcf_readers = [MockVcfReader("PA.vcf"), MockVcfReader("PB.vcf")]
        buffered_readers = merge._create_buffered_readers(vcf_readers)

        self.assertEquals(2, len(buffered_readers))

    def test_create_buffered_readers_modifiesRecords(self):
        file_contents1 = ["##metaheader1\n",
                         '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n',
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|AF:DP|0.24:56|0.01:24")]
        mock_reader1 = MockFileReader("fileA.txt", file_contents1)
        vcf_reader1 = merge.MergeVcfReader(mock_reader1)

        file_contents2 = ["##metaheader1\n",
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\n',
                         self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleNormal|SampleTumor\n"),
                         self.entab("chr2|1|.|A|C|.|.|INFO|DP|32|78")]
        mock_reader2 = MockFileReader("fileB.txt", file_contents2)
        vcf_reader2 = merge.MergeVcfReader(mock_reader2)

        vcf_reader1.store_format_tags("DP", "JX1_DP")
        vcf_reader1.store_format_tags("AF", "AF")
        vcf_reader2.store_format_tags("DP", "JX2_DP")

        buffered_readers = merge._create_buffered_readers([vcf_reader1, vcf_reader2])
        normal_dict = OrderedDict(sorted({"AF":"0.24", "JX1_DP":"56"}.items()))
        tumor_dict = OrderedDict(sorted({"AF":"0.01", "JX1_DP":"24"}.items()))
        expected = OrderedDict(sorted({"fileA|SampleNormal": normal_dict,
                                       "fileA|SampleTumor": tumor_dict}.items()))
        self.assertEquals(expected, buffered_readers[0]._current_element.sample_tag_values)

        normal_dict = OrderedDict(sorted({"JX2_DP":"32"}.items()))
        tumor_dict = OrderedDict(sorted({"JX2_DP":"78"}.items()))
        expected = OrderedDict(sorted({"fileB|SampleNormal": normal_dict,
                                       "fileB|SampleTumor": tumor_dict}.items()))
        self.assertEquals(expected, buffered_readers[1]._current_element.sample_tag_values)

    def test_build_info_tags_sorts(self):
        records = [VcfRecord("1", "42", "A", "C", info="foo"),
                   VcfRecord("1", "43", "A", "C", info="bar")]
        actual_tags = merge._build_info_tags(records)
        self.assertEquals(["bar", "foo"], actual_tags)

    def test_build_info_tags_noDuplicates(self):
        records = [VcfRecord("1", "42", "A", "C", info="foo"),
                   VcfRecord("1", "43", "A", "C", info="bar;foo")]
        actual_tags = merge._build_info_tags(records)
        self.assertEquals(["bar", "foo"], actual_tags)

    def test_build_info_tags_empty(self):
        records = [VcfRecord("1", "42", "A", "C", info="")]
        actual_tags = merge._build_info_tags(records)
        self.assertEquals([], actual_tags)

    def test_build_info_tags_null(self):
        records = [VcfRecord("1", "42", "A", "C", info=".")]
        actual_tags = merge._build_info_tags(records)
        self.assertEquals([], actual_tags)

    def test_build_format_tags_includesModifiedFormatTags(self):
        meta_headers1 = ['##FORMAT=<ID=JX1_DP,Number=1,Type=Integer,Description="Read Depth">',
                        '##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">']
        vcf_reader1 = MockVcfReader(metaheaders=meta_headers1)
        meta_headers2 = ['##FORMAT=<ID=JX2_DP,Number=1,Type=Integer,Description="Read Depth">']
        vcf_reader2 = MockVcfReader(metaheaders=meta_headers2)

        vcf_reader1.store_format_tags("DP", "JX1_DP")
        vcf_reader1.store_format_tags("FA", "FA")
        vcf_reader2.store_format_tags("DP", "JX2_DP")

        vcf_readers = [vcf_reader1, vcf_reader2]
        actual_format_tags = merge._build_format_tags(["^DP$", "^FA$"], vcf_readers)
        expected_format_tags = ["FA", "JX1_DP", "JX2_DP"]

        self.assertEquals(expected_format_tags, actual_format_tags)

    def test_build_format_tags_sameDescriptionOkay(self):
        meta_headers = ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                        '##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">']
        vcf_reader1 = MockVcfReader(metaheaders=meta_headers)
        vcf_reader2 = MockVcfReader(metaheaders=['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">'])

        vcf_reader1.store_format_tags("DP", "DP")
        vcf_reader1.store_format_tags("FA", "FA")
        vcf_reader2.store_format_tags("DP", "DP")

        actual_format_tags = merge._build_format_tags([".*"], [vcf_reader1, vcf_reader2])
        expected_format_tags = ["DP", "FA"]

        self.assertEquals(expected_format_tags, actual_format_tags)

    def test_build_format_tags_filtersAndSorts(self):
        meta_headers = ['##FORMAT=<ID=JQ1>',
                        '##FORMAT=<ID=JQ2>',
                        '##FORMAT=<ID=JQ3>']
        vcf_reader = MockVcfReader(metaheaders=meta_headers)
        vcf_reader.store_format_tags("JQ1", "JQ1")
        vcf_reader.store_format_tags("JQ2", "JQ2")
        vcf_reader.store_format_tags("JQ3", "JQ3")

        format_tags = merge._build_format_tags(["JQ[12]"], [vcf_reader])

        self.assertEquals(["JQ1", "JQ2"], format_tags)

    def test_build_format_tags_uniqueSetAcrossMultipleReaders(self):
        reader1 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1>'])
        reader2 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ2>', '##FORMAT=<ID=JQ3>'])
        reader3 = MockVcfReader(metaheaders=['##INFO=<ID=JQ4>'])

        reader1.store_format_tags("JQ1", "JQ1")
        reader2.store_format_tags("JQ2", "JQ2")
        reader2.store_format_tags("JQ3", "JQ3")
        reader3.store_format_tags("JQ4", "JQ4")

        format_tags = merge._build_format_tags(["JQ[123]"], [reader1, reader2, reader3])

        self.assertEquals(["JQ1", "JQ2", "JQ3"], format_tags)

    def test_build_format_tags_errorIfExcludesAllTags(self):
        reader1 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1>'])
        reader2 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ2>'])

        self.assertRaisesRegexp(utils.JQException,
                                r"The specified format tag regex \[.*\] would exclude all format tags. Review inputs/usage and try again",
                                merge._build_format_tags,
                                ["foo"],
                                [reader1, reader2])

    def test_build_format_tags_warnIfRegexNotFoundInFile(self):
        reader1 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1>'])
        reader2 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ2>'])

        reader1.store_format_tags("JQ1", "JQ1")
        reader2.store_format_tags("JQ2", "JQ2")

        dummy = merge._build_format_tags(["JQ[123]", "foo"], [reader1, reader2])
        actual_log_warnings = test.utils.mock_logger.messages["WARNING"]
        #pylint: disable=anomalous-backslash-in-string
        expected_log_warnings = "In the specified list of regexes \[.*\], the regex \[.*\] does not match any format tags; this expression may be irrelevant."
        self.assertRegexpMatches(actual_log_warnings[0], expected_log_warnings)

    def test_build_contigs(self):
        records = [VcfRecord("1", "42", "A", "C"),
                   VcfRecord("4", "42", "A", "C"),
                   VcfRecord("4", "44", "A", "C"),
                   VcfRecord("15", "42", "A", "C")]
        actual_contigs = merge._build_contigs(records)

        expected_contigs = ["1", "4", "15"]
        self.assertEquals(expected_contigs, actual_contigs)

    def test_compile_metaheaders_preservesExistingMetaheaders(self):
        reader1 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1>'])
        vcf_readers = [reader1]
        all_sample_names = []
        contigs_to_keep = []
        format_tags_to_keep = []
        info_tags_to_keep = []
        existing_metaheaders = ["##existing1", "##existing2"]
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)
        self.assertEquals(3, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertEquals("##existing1", next(metaheaders_iter))
        self.assertEquals("##existing2", next(metaheaders_iter))
        self.assertRegexpMatches(next(metaheaders_iter), "^#CHROM")

    def test_compile_metaheaders_retainsFormatMetaheaders(self):
        reader1 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1,Description="JQA">',
                                             '##FORMAT=<ID=JQ2,Description="JQB">',
                                             '##FORMAT=<ID=JQ3,Description="JQC">'])
        vcf_readers = [reader1]
        all_sample_names = []
        contigs_to_keep = []
        format_tags_to_keep = ["JQ1", "JQ3"]
        info_tags_to_keep = []
        existing_metaheaders = []
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)
        self.assertEquals(3, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(next(metaheaders_iter), "##FORMAT.*JQA")
        self.assertRegexpMatches(next(metaheaders_iter), "##FORMAT.*JQC")
        self.assertRegexpMatches(next(metaheaders_iter), "^#CHROM")

    def test_compile_metaheaders_retainsInfoMetaheaders(self):
        reader1 = MockVcfReader(metaheaders=['##INFO=<ID=JQ1,Description="JQA">',
                                             '##INFO=<ID=JQ2,Description="JQB">',
                                             '##INFO=<ID=JQ3,Description="JQC">'])
        vcf_readers = [reader1]
        all_sample_names = []
        contigs_to_keep = []
        format_tags_to_keep = []
        info_tags_to_keep = ["JQ1", "JQ3", merge._MULT_ALT_TAG]
        existing_metaheaders = []
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)
        self.assertEquals(4, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(next(metaheaders_iter), "##INFO.*JQA")
        self.assertRegexpMatches(next(metaheaders_iter), "##INFO.*JQC")
        self.assertRegexpMatches(next(metaheaders_iter), merge._MULT_ALT_HEADER)
        self.assertRegexpMatches(next(metaheaders_iter), "^#CHROM")

    def test_compile_metaheaders_retainsContigMetaheaders(self):
        reader1 = MockVcfReader(metaheaders=['##contig=<ID=chr1,Description="chromosome 1">',
                                             '##contig=<ID=chr2,Description="chromosome 2">',
                                             '##contig=<ID=chr4,Description="chromosome 2">',
                                             '##contig=<ID=chr11,Description="chromosome 11">'])
        vcf_readers = [reader1]
        all_sample_names = []
        contigs_to_keep = ["chr1", "chr2", "chr11"]
        format_tags_to_keep = []
        info_tags_to_keep = []
        existing_metaheaders = []
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)

        self.assertEquals(4, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(next(metaheaders_iter), "##contig.*chromosome 1")
        self.assertRegexpMatches(next(metaheaders_iter), "##contig.*chromosome 2")
        self.assertRegexpMatches(next(metaheaders_iter), "##contig.*chromosome 11")
        self.assertRegexpMatches(next(metaheaders_iter), "^#CHROM")

    def test_compile_metaheaders_noContigMetaheaders(self):
        reader1 = MockVcfReader(metaheaders=['##foo'])
        vcf_readers = [reader1]
        all_sample_names = []
        contigs_to_keep = ["chr1", "chr2", "chr11"]
        format_tags_to_keep = []
        info_tags_to_keep = []
        existing_metaheaders = []
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)

        self.assertEquals(1, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(next(metaheaders_iter), "^#CHROM")

    def test_compile_metaheaders_addsSampleNamesToColumnHeader(self):
        reader1 = MockVcfReader()
        vcf_readers = [reader1]
        all_sample_names = ["SA", "SB"]
        contigs_to_keep = []
        format_tags_to_keep = []
        info_tags_to_keep = []
        existing_metaheaders = []
        actual_metaheaders = merge._compile_metaheaders(existing_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)
        self.assertEquals(1, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(next(metaheaders_iter), "#CHROM.*SA\tSB")

    def test_compile_metaheaders_ordersHeaders(self):
        reader1 = MockVcfReader(metaheaders=['##FORMAT=<ID=JQ1,Description="JQA">'])
        reader2 = MockVcfReader(metaheaders=['##INFO=<ID=JQ2,Description="JQB">'])

        vcf_readers = [reader1, reader2]
        all_sample_names = []
        contigs_to_keep = []
        format_tags_to_keep = ["JQ1"]
        info_tags_to_keep = ["JQ2"]
        incoming_metaheaders = ["##foo"]
        actual_metaheaders = merge._compile_metaheaders(incoming_metaheaders,
                                                        vcf_readers,
                                                        all_sample_names,
                                                        contigs_to_keep,
                                                        format_tags_to_keep,
                                                        info_tags_to_keep)
        self.assertEquals(4, len(actual_metaheaders))
        metaheaders_iter = iter(actual_metaheaders)
        self.assertRegexpMatches(next(metaheaders_iter), "##foo")
        self.assertRegexpMatches(next(metaheaders_iter), "##INFO.*JQB")
        self.assertRegexpMatches(next(metaheaders_iter), "##FORMAT.*JQA")
        self.assertRegexpMatches(next(metaheaders_iter), "#CHROM.*")

    def test_execute(self):
        vcf_content1 = ('''##source=strelka
##contig=<ID=chr1,Number=1>
##contig=<ID=chr2,Number=1>
##FORMAT=<ID=JQ_Foo1,Number=1,Type=Float,Description="foo">
##FORMAT=<ID=JQ_Bar1,Number=1,Type=Float,Description="bar">
##file1
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleA|SampleB
chr1|1|.|A|C|.|.|INFO|JQ_Foo1:JQ_Bar1|A_1_1:A_1_2|B_1_1:B_1_2
chr1|1|.|A|T|.|.|INFO|JQ_Foo1|A_2|B_2
chr2|1|.|A|C|.|.|INFO|JQ_Foo1:JQ_Bar1|A_3_1:A_3_2|B_3_1:B_3_2
''').replace('|', "\t")
        vcf_content2 = ('''##source=strelka
##contig=<ID=chr1,Number=1>
##contig=<ID=chr2,Number=1>
##file2
##FORMAT=<ID=JQ_Foo2,Number=1,Type=Float,Description="foo">
##FORMAT=<ID=JQ_Bar2,Number=1,Type=Float,Description="bar">
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleA|SampleB
chr1|10|.|A|C|.|.|INFO|JQ_Foo2|C_1_1|D_1_2
chr2|10|.|A|C|.|.|INFO|JQ_Bar2|C_2|D_2
''').replace('|', "\t")

        with TempDirectory() as input_file, TempDirectory() as output_file:
            input_file.write("P1.fileA.vcf", vcf_content1.encode("utf8"))
            input_file.write("P1.fileB.vcf", vcf_content2.encode("utf8"))
            args = Namespace(input=input_file.path,
                             output=os.path.join(output_file.path, "fileC.vcf"),
                             tags=None,
                             include_all=False,
                             include_cells="all",
                             include_rows="all")

            merge.execute(args, ["##execution_header1", "##execution_header2"])

            output_file.check("fileC.vcf")
            with open(os.path.join(output_file.path, "fileC.vcf")) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()

        self.assertEquals(18, len(actual_output_lines))
        actual_lines_iter = iter(actual_output_lines)
        self.assertEquals("##fileformat=VCFv4.1\n", next(actual_lines_iter))
        self.assertEquals("##execution_header1\n", next(actual_lines_iter))
        self.assertEquals("##execution_header2\n", next(actual_lines_iter))
        self.assertRegexpMatches(next(actual_lines_iter), "##jacquard.merge.sample=<Column=1.*>\n")
        self.assertRegexpMatches(next(actual_lines_iter), "##jacquard.merge.sample=<Column=2.*>\n")
        self.assertRegexpMatches(next(actual_lines_iter), "##contig=<ID=chr1,Number=1>\n")
        self.assertRegexpMatches(next(actual_lines_iter), "##contig=<ID=chr2,Number=1>\n")
        self.assertRegexpMatches(next(actual_lines_iter), '##INFO=<ID=JQ_MULT_ALT_LOCUS.*>\n')
        self.assertRegexpMatches(next(actual_lines_iter), '##FORMAT=<ID=JQ_Bar1.*>\n')
        self.assertRegexpMatches(next(actual_lines_iter), '##FORMAT=<ID=JQ_Bar2.*>\n')
        self.assertRegexpMatches(next(actual_lines_iter), '##FORMAT=<ID=JQ_Foo1.*>\n')
        self.assertRegexpMatches(next(actual_lines_iter), '##FORMAT=<ID=JQ_Foo2.*>\n')
        self.assertRegexpMatches(next(actual_lines_iter), "#CHROM\t.*P1|SampleA\tP1|SampleB\n")
        self.assertEquals("chr1\t1\t.\tA\tC\t.\t.\tJQ_MULT_ALT_LOCUS\tJQ_Bar1:JQ_Foo1\tA_1_2:A_1_1\tB_1_2:B_1_1\n", next(actual_lines_iter))
        self.assertEquals("chr1\t1\t.\tA\tT\t.\t.\tJQ_MULT_ALT_LOCUS\tJQ_Foo1\tA_2\tB_2\n", next(actual_lines_iter))
        self.assertEquals("chr1\t10\t.\tA\tC\t.\t.\t.\tJQ_Foo2\tC_1_1\tD_1_2\n", next(actual_lines_iter))
        self.assertEquals("chr2\t1\t.\tA\tC\t.\t.\t.\tJQ_Bar1:JQ_Foo1\tA_3_2:A_3_1\tB_3_2:B_3_1\n", next(actual_lines_iter))
        self.assertEquals("chr2\t10\t.\tA\tC\t.\t.\t.\tJQ_Bar2\tC_2\tD_2\n", next(actual_lines_iter))

    def test_execute_includeFormatIds(self):
        vcf_content1 = ('''##source=strelka
##contig=<ID=chr1,Number=1>
##contig=<ID=chr2,Number=1>
##FORMAT=<ID=JQ_Foo,Number=1,Type=Float,Description="foo">
##FORMAT=<ID=JQ_Bar1,Number=1,Type=Float,Description="bar">
##file1
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SampleA|SampleB
chr1|1|.|A|C|.|.|INFO|JQ_Foo:JQ_Bar1|A_1_1:A_1_2|B_1_1:B_1_2
chr1|1|.|A|T|.|.|INFO|JQ_Foo|A_2|B_2
chr2|1|.|A|C|.|.|INFO|JQ_Foo:JQ_Bar1|A_3_1:A_3_2|B_3_1:B_3_2
''').replace('|', "\t")
        with TempDirectory() as input_file, TempDirectory() as output_file:
            input_file.write("P1.fileA.vcf", vcf_content1.encode("utf8"))
            args = Namespace(input=input_file.path,
                             output=os.path.join(output_file.path, "fileB.vcf"),
                             tags="JQ_.*",
                             include_all=False,
                             include_cells="all",
                             include_rows="all")

            merge.execute(args, ["##extra_header1", "##extra_header2"])

            output_file.check("fileB.vcf")
            with open(os.path.join(output_file.path, "fileB.vcf")) as actual_output_file:
                actual_output_lines = actual_output_file.readlines()

        expected_output_headers = ["##fileformat=VCFv4.1\n",
                                   "##extra_header1\n",
                                   "##extra_header2\n",
                                   "##jacquard.merge.sample=<Column=1,Name=P1|SampleA,Source=P1.fileA.vcf>\n",
                                   "##jacquard.merge.sample=<Column=2,Name=P1|SampleB,Source=P1.fileA.vcf>\n",
                                   '##contig=<ID=chr1,Number=1>\n',
                                   '##contig=<ID=chr2,Number=1>\n',
                                   '##INFO=<ID=JQ_MULT_ALT_LOCUS,Number=0,Type=Flag,Description="More than one alt allele was seen at this locus.">\n',
                                   '##FORMAT=<ID=JQ_Bar1,Number=1,Type=Float,Description="bar">\n',
                                   '##FORMAT=<ID=JQ_Foo,Number=1,Type=Float,Description="foo">\n',
                                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1|SampleA\tP1|SampleB\n"]

        self.assertEquals(expected_output_headers, actual_output_lines[0:len(expected_output_headers)])

    def test_sort_readers_orderedVcfsPassThrough(self):
        record1 = vcf.VcfRecord("chr1", "42", "A", "C")
        record2 = vcf.VcfRecord("chr2", "42", "A", "C")
        record3 = vcf.VcfRecord("chr3", "42", "A", "C")
        vcf_readerA = MockVcfReader(records=[record1, record2, record3])
        vcf_readerB = MockVcfReader(records=[record1, record2, record3])

        input_readers = [vcf_readerA, vcf_readerB]
        with TempDirectory() as temp_dir:
            actual_readers = merge._sort_readers(list(input_readers),
                                                 temp_dir.path)

            self.assertEquals(actual_readers, input_readers)

    def test_sort_readers_unsortedReturnsMergeVcfReader(self):
        record1 = vcf.VcfRecord("chr1", "42", "A", "C")
        record2 = vcf.VcfRecord("chr3", "42", "A", "C")
        record3 = vcf.VcfRecord("chr2", "42", "A", "C")
        vcf_readerA = MockVcfReader(records=[record1, record2, record3])
        vcf_readerB = MockVcfReader(records=[record1, record2, record3])

        input_readers = [vcf_readerA, vcf_readerB]
        with TempDirectory() as temp_dir:
            actual_readers = merge._sort_readers(list(input_readers),
                                                 temp_dir.path)

            self.assertEquals("MergeVcfReader", type(actual_readers[0]).__name__)
            self.assertEquals("MergeVcfReader", type(actual_readers[1]).__name__)

    def test_sort_readers_vcfsResortedAsNecessary(self):
        #pylint: disable=too-many-locals
        record1 = vcf.VcfRecord("chr1", "42", "A", "C")
        record2 = vcf.VcfRecord("chr2", "42", "A", "C")
        record3 = vcf.VcfRecord("chr3", "42", "A", "C")
        vcf_readerA = MockVcfReader(records=[record1, record2, record3])
        input_metaheaders = ["##foo", "##bar"]
        input_column_header = self.entab("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SX|SY")
        vcf_readerB = MockVcfReader(input_filepath="unsorted.vcf",
                                    metaheaders=input_metaheaders,
                                    column_header=input_column_header,
                                    records=[record3, record1, record2])

        input_readers = [vcf_readerA, vcf_readerB]
        with TempDirectory() as temp_dir:
            temp_dir.write("foo.vcf", b"foo\nbar")
            output_filepath = os.path.join(temp_dir.path, "foo.vcf")
            actual_readers = merge._sort_readers(list(input_readers),
                                                 output_filepath)

            self.assertNotEquals(actual_readers, input_readers)
            self.assertIn(vcf_readerA, actual_readers)
            self.assertNotIn(vcf_readerB, actual_readers)
            self.assertEquals(2, len(actual_readers))
            actual_readerB = actual_readers[1]
            actual_readerB.open()
            actual_metaheaders = list(actual_readerB.metaheaders)
            actual_column_header = actual_readerB.column_header
            actual_records = [vcf_record for vcf_record in actual_readerB.vcf_records()]
            actual_readerB.close()


        self.assertEquals(input_metaheaders, actual_metaheaders)
        self.assertEquals(input_column_header, actual_column_header)
        self.assertEquals(record1.text(), actual_records[0].text())
        self.assertEquals(record2.text(), actual_records[1].text())
        self.assertEquals(record3.text(), actual_records[2].text())
        actual_log_infos = test.utils.mock_logger.messages["INFO"]
        self.assertEquals(3, len(actual_log_infos))
        self.assertRegexpMatches(actual_log_infos[0], r"Checking sort order of \[vcfName\] \(1/2\)")
        self.assertRegexpMatches(actual_log_infos[1], r"Checking sort order of \[unsorted.vcf\] \(2/2\)")
        self.assertRegexpMatches(actual_log_infos[2], r"Sorting vcf \[unsorted.vcf\] \(1/1\)")
        actual_log_debugs = test.utils.mock_logger.messages["DEBUG"]
        self.assertEquals(1, len(actual_log_debugs))
        self.assertRegexpMatches(actual_log_debugs[0],
                                 r"VCF file:chrom:pos \[unsorted.vcf:chr1:42\] is out of order")


class MergeFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_merge(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))

            module_testdir = os.path.join(test_dir,
                                          "functional_tests",
                                          "02_merge")
            input_file = os.path.join(module_testdir, "input")
            output_file = os.path.join(output_dir.path, "merged.vcf")

            command = ["merge", input_file, output_file, "--force"]
            expected_file = os.path.join(module_testdir,
                                        "benchmark",
                                        "merged.vcf")

            self.assertCommand(command, expected_file)

    def test_merge_unsorted(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))

            module_testdir = os.path.join(test_dir,
                                          "functional_tests",
                                          "02_merge_unsorted")
            input_file = os.path.join(module_testdir, "input")
            output_file = os.path.join(output_dir.path, "merged.vcf")

            command = ["merge", input_file, output_file, "--force"]
            expected_file = os.path.join(module_testdir,
                                        "benchmark",
                                        "merged.vcf")

            self.assertCommand(command, expected_file)

class BufferedReaderTestCase(test_case.JacquardBaseTestCase):
    def test_get_sample_info_advancesCurrentElementWhenMatched(self):
        reader = [1, 5, 10, 15]
        buffered_reader = merge._BufferedReader(iter(reader))

        self.assertEquals(None, buffered_reader.next_if_equals(0))
        self.assertEquals(None, buffered_reader.next_if_equals(5))
        self.assertEquals(1, buffered_reader.next_if_equals(1))
        self.assertEquals(None, buffered_reader.next_if_equals(1))
        self.assertEquals(5, buffered_reader.next_if_equals(5))
        self.assertEquals(10, buffered_reader.next_if_equals(10))
        self.assertEquals(15, buffered_reader.next_if_equals(15))
        self.assertEquals(None, buffered_reader.next_if_equals(15))
        self.assertEquals(None, buffered_reader.next_if_equals(42))

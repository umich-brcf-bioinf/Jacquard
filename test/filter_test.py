# pylint: disable=line-too-long, invalid-name, global-statement, star-args, too-many-public-methods, too-few-public-methods
from __future__ import absolute_import

from argparse import Namespace
import os

from testfixtures import TempDirectory

import jacquard.commands.filter as filter
import jacquard.logger
import test.mock_logger
import test.test_case as test_case
from test.vcf_test import MockWriter, MockVcfReader, MockVcfRecord


#TODO: (cgates): These tests should start using mocked readers/writers and stop using the file system
VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\tsampleB\n"

class FilterSomaticTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(FilterSomaticTestCase, self).setUp()
        filter.logger = test.mock_logger

    def tearDown(self):
        filter.logger = jacquard.logger
        super(FilterSomaticTestCase, self).tearDown()

    def test_validate_arguments(self):
        with TempDirectory() as input_file, TempDirectory() as output_file:
            input_file.write("A.normalized.vcf", b"##source=strelka\n#colHeader")
            input_file.write("B.normalized.vcf", b"##source=strelka\n#colHeader")
            args = Namespace(input=input_file.path,
                             output=output_file.path)

            actual_readers_to_writers, out_files = filter._validate_arguments(args)
            self.assertEquals(2, len(actual_readers_to_writers))
            self.assertEquals(0, len(out_files))
            self.assertRegexpMatches(list(actual_readers_to_writers.values())[0].output_filepath, "HCsomatic.vcf")

    def test_predict_output(self):
        with TempDirectory() as input_file:
            input_file.write("A.normalized.jacquardTags.vcf", b"##source=strelka\n#colHeader")
            input_file.write("B.normalized.jacquardTags.vcf", b"##source=strelka\n#colHeader")
            args = Namespace(input=input_file.path)

            desired_output_files = filter._predict_output(args)
            expected_desired_output_files = set(["A.normalized.jacquardTags.HCsomatic.vcf",
                                                 "B.normalized.jacquardTags.HCsomatic.vcf"])

            self.assertEquals(expected_desired_output_files, desired_output_files)

            self.assertEquals(expected_desired_output_files, desired_output_files)

    def test_find_positions(self):
        vcf_readers = [MockVcfReader(input_filepath="file1",
                                     records=[MockVcfRecord("1", "2352", "A", "G", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])]),
                       MockVcfReader(input_filepath="file2",
                                     records=[MockVcfRecord("1", "2353", "A", "GT", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])])]
        somatic_positions = filter._find_positions(vcf_readers, None, None)
        self.assertEquals({'1^2353^A': 1, '1^2352^A': 1}, somatic_positions)

    def test_find_positions_distinctRefs(self):
        vcf_readers = [MockVcfReader(input_filepath="file1",
                                     records=[MockVcfRecord("1", "2352", "AT", "G", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])]),
                       MockVcfReader(input_filepath="file2",
                                     records=[MockVcfRecord("1", "2352", "A", "G", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])])]
        somatic_positions = filter._find_positions(vcf_readers, None, None)
        self.assertEquals({'1^2352^A': 1, '1^2352^AT': 1}, somatic_positions)

    def test_find_positions_filterJQExclude(self):
        vcf_readers = [MockVcfReader(input_filepath="file1",
                                     records=[MockVcfRecord("1", "2352", "A", "GT", vcf_filter="JQ_EXCLUDE", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])]),
                       MockVcfReader(input_filepath="file2",
                                     records=[MockVcfRecord("1", "2353", "A", "G", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])])]
        somatic_positions = filter._find_positions(vcf_readers, None, None)
        self.assertEquals({'1^2353^A': 1}, somatic_positions)

    def test_find_positions_filterJQExclude_messages(self):
        vcf_readers = [MockVcfReader(input_filepath="file1",
                                     records=[MockVcfRecord("1", "2352", "A", "GT", vcf_filter="JQ_EXCLUDE", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])]),
                       MockVcfReader(input_filepath="file2",
                                     records=[MockVcfRecord("1", "2353", "A", "GC", vcf_filter="JQ_EXCLUDE", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])]),
                       MockVcfReader(input_filepath="file3",
                                     records=[MockVcfRecord("1", "2354", "A", "G", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])])]
        somatic_positions = filter._find_positions(vcf_readers, None, None)
        self.assertEquals({'1^2354^A': 1}, somatic_positions)

        actual_log_debugs = test.mock_logger.messages["DEBUG"]
        self.assertIn('Removed [2] problematic mockCaller variant records with filter=JQ_EXCLUDE', actual_log_debugs)
        actual_log_warnings = test.mock_logger.messages["WARNING"]
        self.assertIn("A total of [2] problematic variant records failed Jacquard's filters. See output and log for details.", actual_log_warnings)

    def test_find_positions_noSomaticPositions(self):
        vcf_readers = [MockVcfReader(input_filepath="file1",
                                     records=[MockVcfRecord("1", "2352", "A", "GT", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:0", "52:0"])]),
                       MockVcfReader(input_filepath="file2",
                                     records=[MockVcfRecord("1", "2353", "A", "GC", vcf_format="DP", samples=["234", "52"])])]
        dummy = filter._find_positions(vcf_readers, None, None)
        actual_log_warnings = test.mock_logger.messages["WARNING"]
        self.assertIn('Input file [file1] has no high-confidence somatic variants.', actual_log_warnings)
        self.assertIn('Input file [file2] has no high-confidence somatic variants.', actual_log_warnings)
        self.assertIn('[2] VCF file(s) had no high-confidence somatic variants. See log for details.', actual_log_warnings)

    def test_find_positions_includePassedVariants(self):
        vcf_readers = [MockVcfReader(input_filepath="file1",
                                     records=[MockVcfRecord("1", "2352", "A", "G", vcf_filter="PASS", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])]),
                       MockVcfReader(input_filepath="file2",
                                     records=[MockVcfRecord("1", "2353", "A", "GT", vcf_filter="FAIL", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])])]
        somatic_positions = filter._find_positions(vcf_readers, "passed", None)
        self.assertEquals({'1^2352^A': 1}, somatic_positions)

    def test_find_positions_includeSomaticVariants(self):
        vcf_readers = [MockVcfReader(input_filepath="file1",
                                     records=[MockVcfRecord("1", "2352", "A", "G", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])]),
                       MockVcfReader(input_filepath="file2",
                                     records=[MockVcfRecord("1", "2353", "A", "GT", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:0", "52:0"])])]
        somatic_positions = filter._find_positions(vcf_readers, "somatic", None)
        self.assertEquals({'1^2352^A': 1}, somatic_positions)

    def test_find_positions_includePassedLoci(self):
        vcf_readers = [MockVcfReader(input_filepath="file1",
                                     records=[MockVcfRecord("1", "2352", "A", "G", vcf_filter="PASS", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])]),
                       MockVcfReader(input_filepath="file2",
                                     records=[MockVcfRecord("1", "2353", "A", "GT", vcf_filter="FAIL", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])])]
        somatic_positions = filter._find_positions(vcf_readers, None, "all_passed")
        self.assertEquals({'1^2352^A': 1}, somatic_positions)

    def test_find_positions_includeSomaticLoci(self):
        vcf_readers = [MockVcfReader(input_filepath="file1",
                                     records=[MockVcfRecord("1", "2352", "A", "G", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:1", "52:1"])]),
                       MockVcfReader(input_filepath="file2",
                                     records=[MockVcfRecord("1", "2353", "A", "GT", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:0", "52:0"])])]
        somatic_positions = filter._find_positions(vcf_readers, None, "all_somatic")
        self.assertEquals({'1^2352^A': 1}, somatic_positions)

    def test_writePositions_outputFileForEachInputFile(self):
        vcf_reader1 = MockVcfReader(input_filepath="mutect.vcf",
                                    records=[MockVcfRecord("1", "32", "A", "G", vcf_format="DP:JQ_HC_SOM_MT", samples=["234:1", "52:1"])])
        vcf_reader2 = MockVcfReader(input_filepath="varscan.vcf",
                                    records=[MockVcfRecord("1", "35", "A", "GT", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:0", "52:0"]),
                                             MockVcfRecord("1", "32", "A", "GT", vcf_format="DP:JQ_HC_SOM_VS", samples=["234:0", "52:0"])])
        vcf_writer1 = MockWriter(output_filepath="mutect.HCsomatic.vcf")
        vcf_writer2 = MockWriter(output_filepath="varscan.HCsomatic.vcf")
        readers_to_writers = {vcf_reader1: vcf_writer1, vcf_reader2: vcf_writer2}

        positions = {"1^32^A"}
        execution_context = ["##foo", "##bar"]

        filter._write_positions(readers_to_writers, positions, execution_context, None)

        actual_log_warnings = test.mock_logger.messages["INFO"]
        self.assertIn("Filtered to [2] calls in high-confidence loci.", actual_log_warnings)

    def test_writeOutput(self):
        mock_writer = MockWriter()
        headers = ["#foo", "#bar"]
        actual_sorted_variants = ["123", "456"]

        filter._write_output(mock_writer, headers, actual_sorted_variants)
        actualLines = mock_writer.lines()

        self.assertEqual("#foo", actualLines[0])
        self.assertEqual("#bar", actualLines[1])
        self.assertEqual("123", actualLines[2])
        self.assertEqual("456", actualLines[3])

class FilterFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_filter(self):
        with TempDirectory() as output_file:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "03_filter")
            input_file = os.path.join(module_testdir, "input")

            command = ["filter", input_file, output_file.path, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

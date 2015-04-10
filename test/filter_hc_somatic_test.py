# pylint: disable=line-too-long, invalid-name, global-statement, star-args, too-many-public-methods, too-few-public-methods
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import os

from testfixtures import TempDirectory

import jacquard.commands.filter_hc_somatic as filter_hc_somatic
import jacquard.logger
import test.mock_logger
import test.test_case as test_case


#TODO: (cgates): These tests should start using mocked readers/writers and stop using the file system
VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\tsampleB\n"

class FilterSomaticTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(FilterSomaticTestCase, self).setUp()
        filter_hc_somatic.logger = test.mock_logger

    def tearDown(self):
        filter_hc_somatic.logger = jacquard.logger
        super(FilterSomaticTestCase, self).tearDown()

    def test_validate_arguments(self):
        with TempDirectory() as input_file, TempDirectory() as output_file:
            input_file.write("A.normalized.vcf", b"##source=strelka\n#colHeader")
            input_file.write("B.normalized.vcf", b"##source=strelka\n#colHeader")
            args = Namespace(input=input_file.path,
                             output=output_file.path)

            actual_readers_to_writers, out_files = filter_hc_somatic._validate_arguments(args)
            self.assertEquals(2, len(actual_readers_to_writers))
            self.assertEquals(0, len(out_files))
            self.assertRegexpMatches(list(actual_readers_to_writers.values())[0].output_filepath,
                                     "HCsomatic.vcf")

    def test_predict_output(self):
        with TempDirectory() as input_file:
            input_file.write("A.normalized.jacquardTags.vcf", b"##source=strelka\n#colHeader")
            input_file.write("B.normalized.jacquardTags.vcf", b"##source=strelka\n#colHeader")
            args = Namespace(input=input_file.path)

            desired_output_files = filter_hc_somatic._predict_output(args)
            expected_desired_output_files = set(["A.normalized.jacquardTags.HCsomatic.vcf",
                                                 "B.normalized.jacquardTags.HCsomatic.vcf"])

            self.assertEquals(expected_desired_output_files, desired_output_files)

            self.assertEquals(expected_desired_output_files, desired_output_files)

    def test_find_somatic_positions(self):
        with TempDirectory() as input_file:
            input_file.write("A.snp.vcf",
                             (b"##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
                              b"INFO\tFORMAT\tNORMAL\tTUMOR\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP\t234\n1\t2352"
                              b"\t.\tA\tG\t.\t.\tfoo\tDP:JQ_HC_SOM_VS\t234:1\t52:1\n"))
            input_file.write("A.indel.vcf",
                             (b"##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                              b"\tFORMAT\tNORMAL\tTUMOR\n1\t2353\t.\tA\tGT\t.\t.\tfoo\tDP:JQ_HC_SOM_VS\t234:1\t52:1\n"))

            file1 = os.path.join(input_file.path, "A.snp.vcf")
            file2 = os.path.join(input_file.path, "A.indel.vcf")

            somatic_positions, somatic_positions_header = filter_hc_somatic._find_somatic_positions([file1, file2])
            self.assertEqual({'1^2353^A': 1, '1^2352^A': 1}, somatic_positions)
            self.assertEqual("##jacquard.filterHCSomatic.total_highConfidence_somatic_positions=2\n", somatic_positions_header)

    def test_find_somatic_positions_distinctRefs(self):
        with TempDirectory() as input_file:
            fileA_contents = ("##source=VarScan2\n"
                              "#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR\n"
                              "1|10|.|A|G|.|.foo|DP|234\n"
                              "1|20|.|AT|G|.|.|foo|DP:JQ_HC_SOM_VS|234:1|52:1\n")
            fileB_contents = ("##source=VarScan2\n"
                              "#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR\n"
                              "1|20|.|A|GT|.|.|foo|DP:JQ_HC_SOM_VS|234:1|52:1\n")
            input_file.write("A.vcf", self.entab(fileA_contents).encode("utf8"))
            input_file.write("B.vcf", self.entab(fileB_contents).encode("utf8"))

            input_fileA = os.path.join(input_file.path, "A.vcf")
            input_fileB = os.path.join(input_file.path, "B.vcf")
            somatic_positions, somatic_positions_header = filter_hc_somatic._find_somatic_positions([input_fileA, input_fileB])

            self.assertEqual({'1^20^AT': 1, '1^20^A': 1}, somatic_positions)
            self.assertEqual("##jacquard.filterHCSomatic.total_highConfidence_somatic_positions=2\n", somatic_positions_header)

    def test_filterJQExclude(self):
        with TempDirectory() as input_file:
            input_file.write("A.snp.vcf",
                             (b"##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
                              b"INFO\tFORMAT\tNORMAL\tTUMOR\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP\t234\n1\t2352"+
                              b"\t.\tA\tG\t.\tJQ_EXCLUDE\tfoo\tDP:JQ_HC_SOM_VS\t234:1\t52:1\n"))
            input_file.write("A.indel.vcf",
                             (b"##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                              b"\tFORMAT\tNORMAL\tTUMOR\n1\t2353\t.\tA\tGT\t.\t.\tfoo\tDP:JQ_HC_SOM_VS\t234:1\t52:1\n"))

            file1 = os.path.join(input_file.path, "A.snp.vcf")
            file2 = os.path.join(input_file.path, "A.indel.vcf")

            somatic_positions, somatic_positions_header = filter_hc_somatic._find_somatic_positions([file1, file2])
            self.assertEqual({'1^2353^A': 1}, somatic_positions)
            self.assertEqual("##jacquard.filterHCSomatic.total_highConfidence_somatic_positions=1\n", somatic_positions_header)

            input_file.cleanup()

    def test_filterJQExclude_messages(self):
        with TempDirectory() as input_file:
            input_file.write("A.snp.vcf",
                             (b"##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
                              b"INFO\tFORMAT\tNORMAL\tTUMOR\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP\t234\n1\t2352"
                              b"\t.\tA\tG\t.\tJQ_EXCLUDE\tfoo\tDP:JQ_HC_SOM_VS\t234:1\t52:1\n"))
            input_file.write("A.indel.vcf",
                             (b"##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                              b"\tFORMAT\tNORMAL\tTUMOR\n1\t2353\t.\tA\tGT\t.\t.\tfoo\tDP:JQ_HC_SOM_VS\t234:1\t52:1\n"))
            input_file.write("A.snvs.vcf",
                             (b"##source=strelka\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
                              b"INFO\tFORMAT\tSAMPLE\n\t2352\t.\tA\tG\t.\t.\tfoo\tDP\t234\n1\t2352"
                              b"\t.\tA\tG\t.\tJQ_EXCLUDE\tfoo\tDP:JQ_HC_SOM_SK\t234:1"))
            input_file.write("A.indels.vcf",
                             (b"##source=strelka\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                              b"\tFORMAT\tSAMPLE\n1\t2353\t.\tA\tGT\t.\t.\tfoo\tDP:JQ_HC_SOM_SK\t234:1\n"))

            files = [os.path.join(input_file.path, "A.snp.vcf"),
                     os.path.join(input_file.path, "A.indel.vcf"),
                     os.path.join(input_file.path, "A.snvs.vcf"),
                     os.path.join(input_file.path, "A.indels.vcf")]

            filter_hc_somatic._find_somatic_positions(files)

            actual_log_debugs = test.mock_logger.messages["DEBUG"]
            self.assertIn('Removed [1] problematic VarScan variant records with filter=JQ_EXCLUDE', actual_log_debugs)
            self.assertIn('Removed [1] problematic Strelka variant records with filter=JQ_EXCLUDE', actual_log_debugs)
            actual_log_warnings = test.mock_logger.messages["WARNING"]
            self.assertIn("A total of [2] problematic variant records failed Jacquard's filters. See output and log for details.", actual_log_warnings)

    def test_findSomaticPositions_invalidInput(self):
        with TempDirectory() as input_file:
            input_file.write("A.vcf",
                             b"##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP\t234\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP\t234:1\n")
            input_file.write("B.vcf",
                             b"##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n1\t2353\t.\tA\tGT\t.\t.\tfoo\tDP:JQ_SOM_VS\t234:1\n")

            file1 = os.path.join(input_file.path, "A.vcf")
            file2 = os.path.join(input_file.path, "B.vcf")

            filter_hc_somatic._find_somatic_positions([file1, file2])
            actual_log_warnings = test.mock_logger.messages["WARNING"]
            self.assertIn('Input file [A.vcf] has no high-confidence somatic variants.', actual_log_warnings)
            self.assertIn('Input file [B.vcf] has no high-confidence somatic variants.', actual_log_warnings)
            self.assertIn('[2] VCF file(s) had no high-confidence somatic variants. See log for details.', actual_log_warnings)

    def test_writeSomatic_outputFileForEachInputFile(self):
        with TempDirectory() as input_file, TempDirectory() as output_file:
            input1 = input_file.write('mutect.vcf', "##MuTect\n{}1\t32\t.\tA\tT\t.\tPASS\tFOO\tDP\t30".format(VCF_HEADER).encode("utf8"))
            input2 = input_file.write('varscan.vcf', "##source=VarScan2\n{}1\t35\t.\tA\tG\t.\tPASS\tFOO\tDP\t30".format(VCF_HEADER).encode("utf8"))

            in_files = [input1, input2]
            somatic_positions = {"1^32^A"}
            execution_context = ["##foo", "##bar"]
            filter_hc_somatic._write_somatic(in_files, output_file.path, somatic_positions, execution_context)

            self.assertEqual(["mutect.HCsomatic.vcf", "varscan.HCsomatic.vcf"], output_file.actual())
            actual_log_warnings = test.mock_logger.messages["INFO"]
            self.assertIn("Filtered to [1] calls in high-confidence loci.", actual_log_warnings)

    def test_sort_sortHeaders(self):
        headers = ["##foo", "##bar", "#CHROM", "##baz"]
        sorted_headers = filter_hc_somatic._sort_headers(headers)
        expected_sorted_headers = ["##bar", "##baz", "##foo", "#CHROM"]
        self.assertEqual(expected_sorted_headers, sorted_headers)


    def test_writeOutput(self):
        mock_writer = MockWriter()
        headers = ["#foo", "#bar"]
        actual_sorted_variants = ["123", "456"]

        filter_hc_somatic._write_output(mock_writer, headers, actual_sorted_variants)
        actualLines = mock_writer.lines()

        self.assertEqual("#foo", actualLines[0])
        self.assertEqual("#bar", actualLines[1])
        self.assertEqual("123", actualLines[2])
        self.assertEqual("456", actualLines[3])


class MockWriter(object):
    def __init__(self):
        self._content = []
        self.wasClosed = False

    def write(self, content):
        self._content.extend(content.splitlines())

    def lines(self):
        return self._content

    def close(self):
        self.wasClosed = True

class FilterHCSomaticFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_filter_hc_somatic(self):
        with TempDirectory() as output_file:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir, "functional_tests", "02_filter_hc_somatic")
            input_file = os.path.join(module_testdir, "input")

            command = ["filter_hc_somatic", input_file, output_file.path, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

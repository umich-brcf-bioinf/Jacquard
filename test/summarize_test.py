#pylint: disable=invalid-name, global-statement, line-too-long
#pylint: disable=too-many-public-methods, too-few-public-methods, unused-argument
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import os

from testfixtures import TempDirectory

import jacquard.summarize as summarize
import jacquard.utils.vcf as vcf
import test.utils.test_case as test_case
from test.utils.vcf_test import MockFileWriter, MockVcfReader


#TODO (cgates): The module summarize is not adequately unit-tested
class MockSummarizeCaller(object):
    def __init__(self, metaheaders_list=None):
        if metaheaders_list:
            self.metaheaders_list = metaheaders_list
        else:
            self.metaheaders_list = []
        self.add_tags_called = False

    def add_tags(self, vcf_record):
        self.add_tags_called = True

    def get_metaheaders(self):
        return self.metaheaders_list

class SummarizeTestCase(test_case.JacquardBaseTestCase):
    def test_write_metaheaders(self):
        file_writer = MockFileWriter()
        vcf_reader = MockVcfReader(column_header="#header")
        caller = MockSummarizeCaller(["##summarize_metaheader"])
        summarize._write_metaheaders(caller,
                                     vcf_reader,
                                     file_writer,
                                     ["##execution_context"])
        expected = ["##execution_context",
                    "##metaheaders",
                    "##summarize_metaheader",
                    "#header"]
        self.assertEquals(expected, file_writer.lines())

    def test_add_summarize_tags(self):
        file_writer = MockFileWriter()
        vcf_record = vcf.VcfRecord("chr1", "42", "A", "C")
        vcf_reader = MockVcfReader(records=[vcf_record])
        caller = MockSummarizeCaller()

        summarize._add_tags(caller, vcf_reader, file_writer)

        self.assertTrue(caller.add_tags_called)

    def test_execute(self):
        input_data = self.entab(\
"""##blah\n#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SAMPLE
1|42|.|A|G|.|PASS|INFO|JQ_VS_AF:JQ_MT_AF:JQ_VS_DP:JQ_MT_DP|0.2:0.4:30:45""")
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("foo.vcf", input_data.encode("utf8"))
            input_file = os.path.join(input_dir.path, "foo.vcf")
            output_file = os.path.join(output_dir.path, "baz.vcf")
            args = Namespace(input=input_file,
                             output=output_file,
                             column_specification=None)
            summarize.execute(args, ["##foo"])
            output_dir.check("baz.vcf")

class SummarizeFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_summarize(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir,
                                          "functional_tests",
                                          "03_summarize")
            input_dir = os.path.join(module_testdir,
                                     "input",
                                     "tiny_strelka.merged.vcf")
            output_file = os.path.join(output_dir.path, "summarized.vcf")

            command = ["summarize", input_dir, output_file, "--force"]
            expected_file = os.path.join(module_testdir,
                                         "benchmark",
                                         "summarized.vcf")

            self.assertCommand(command, expected_file)


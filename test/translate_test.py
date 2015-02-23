#pylint: disable=line-too-long,too-many-public-methods,invalid-name
#pylint: disable=missing-docstring,protected-access,too-few-public-methods
#pylint: disable=too-many-arguments,too-many-instance-attributes
from argparse import Namespace
import os

from testfixtures import TempDirectory

from jacquard import __version__, vcf
import jacquard.translate as translate
import test.test_case as test_case
from test.vcf_test import MockVcfReader, MockTag, MockWriter

class TranslateTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(TranslateTestCase, self).setUp()

    def test_get_required_input_output_types(self):
        self.assertEquals(("directory", "directory"),
                          translate.get_required_input_output_types())

    def test_report_prediction(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf", "##source=strelka\n#colHeader")
            input_dir.write("B.vcf", "##source=strelka\n#colHeader")
            input_dir.write("B.hpfilter.pass", "##source=strelka\n#colHeader")
            args = Namespace(input=input_dir.path)

            desired_output_files = translate.report_prediction(args)
            expected_desired_output_files = set(["A.translatedTags.vcf",
                                                 "B.translatedTags.vcf"])

            self.assertEquals(expected_desired_output_files, desired_output_files)

    def test_translate_files(self):
        record = vcf.VcfRecord("chr1", "42", "A", "C",
                               sample_tag_values={"SA":{}, "SB":{}})
        reader = MockVcfReader(metaheaders=["##metaheader1",
                                            "##metaheader2"],
                               records=[record],
                               sample_names=["SA", "SB"])
        writer = MockWriter()
        execution_context = []
        new_tags = [MockTag("TAG1", {"SA":42, "SB":43}, metaheader="##newTag1"),
                    MockTag("TAG2", {"SA":420, "SB":430}, metaheader="##newTag2")]
        translate._translate_files(reader,
                                   new_tags,
                                   execution_context,
                                   writer)

        self.assertTrue(reader.opened)
        self.assertTrue(writer.opened)
        expected = ['##metaheader1',
                    '##metaheader2',
                    '##newTag1',
                    '##newTag2',
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR']
        self.assertEquals(expected, writer._content[0:5])
        self.assertRegexpMatches(writer._content[5], "TAG1:TAG2")
        self.assertRegexpMatches(writer._content[5], "42:420")
        self.assertRegexpMatches(writer._content[5], "43:430")

        self.assertTrue(reader.closed)
        self.assertTrue(writer.closed)


    def test_translate_files_sortsRecords(self):
        writer = MockWriter()
        record1 = vcf.VcfRecord("chr1", "42", "A", "C")
        record2 = vcf.VcfRecord("chr2", "42", "A", "C")
        record3 = vcf.VcfRecord("chr3", "42", "A", "C")
        vcf_reader = MockVcfReader(metaheaders=[],
                                   records=[record3, record1, record2])
        execution_context = []
        new_tags = []

        translate._translate_files(vcf_reader,
                                   new_tags,
                                   execution_context,
                                   writer)

        actual_lines= iter(writer.lines())
        self.assertRegexpMatches(actual_lines.next(), "^#CHROM")
        self.assertRegexpMatches(actual_lines.next(), "^chr1")
        self.assertRegexpMatches(actual_lines.next(), "^chr2")
        self.assertRegexpMatches(actual_lines.next(), "^chr3")

    def test_translate_write_metaheaders_addsExecutionMetaheaders(self):
        writer = MockWriter()
        reader = MockVcfReader(metaheaders=["##mockCallerMetaheader1"],
                               column_header="#CHROM\tPOS\tREF\tALT\tStuff")
        execution_context = ["##foo1=bar",
                             "##foo2=baz"]
        new_tags = [MockTag(metaheader="##newTag1"),
                    MockTag(metaheader="##newTag2")]
        translate._write_headers(reader, new_tags, execution_context, writer)
        expected_headers = ["##mockCallerMetaheader1",
                            "##foo1=bar",
                            "##foo2=baz",
                            "##newTag1",
                            "##newTag2",
                            "#CHROM\tPOS\tREF\tALT\tStuff"]
        self.assertEquals(expected_headers, writer._content)


class ExcludeMalformedRefTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEquals('##FILTER=<ID=JQ_EXCLUDE_MALFORMED_REF,Description="The format of the reference value for this variant record does not comply with VCF standard.",Source="Jacquard",Version="{}">'.format(__version__),
                          translate._ExcludeMalformedRef().metaheader)

    def test_add_tag_value_validRefNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "C", vcf_filter="PASS")
        translate._ExcludeMalformedRef().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_validIndelRefNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "ACGT", "C", vcf_filter="PASS")
        translate._ExcludeMalformedRef().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_validIndelRefEdgecaseNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "ACGTNacgtn", "C", vcf_filter="PASS")
        translate._ExcludeMalformedRef().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_invalidRefReplacesFilter(self):
        record = vcf.VcfRecord("chr1", "42", "X", "C", vcf_filter="PASS")
        translate._ExcludeMalformedRef().add_tag_values(record)
        self.assertEquals("JQ_EXCLUDE_MALFORMED_REF", record.filter)

    def test_add_tag_value_invalidIndelReplacesFilter(self):
        record = vcf.VcfRecord("chr1", "42", "XYZ", "C", vcf_filter="PASS")
        translate._ExcludeMalformedRef().add_tag_values(record)
        self.assertEquals("JQ_EXCLUDE_MALFORMED_REF", record.filter)

class ExcludeMalformedAltTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEquals('##FILTER=<ID=JQ_EXCLUDE_MALFORMED_ALT,Description="The the format of the alternate allele value for this variant record does not comply with VCF standard.",Source="Jacquard",Version={}>'.format(__version__),
                          translate._ExcludeMalformedAlt().metaheader)

    def test_add_tag_value_validAltNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "C", vcf_filter="PASS")
        translate._ExcludeMalformedAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_validIndelAltNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "AC,GT", vcf_filter="PASS")
        translate._ExcludeMalformedAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_validIndelAltEdgecaseNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "ACGTNacgtn,", vcf_filter="PASS")
        translate._ExcludeMalformedAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_invalidAltReplacesFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "X", vcf_filter="PASS")
        translate._ExcludeMalformedAlt().add_tag_values(record)
        self.assertEquals("JQ_EXCLUDE_MALFORMED_ALT", record.filter)

    def test_add_tag_value_invalidIndelReplacesFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "XYZ", vcf_filter="PASS")
        translate._ExcludeMalformedAlt().add_tag_values(record)
        self.assertEquals("JQ_EXCLUDE_MALFORMED_ALT", record.filter)

    def test_add_tag_value_missingAltBothReplacesFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", ".*", vcf_filter="PASS")
        translate._ExcludeMalformedAlt().add_tag_values(record)
        self.assertEquals("JQ_EXCLUDE_MALFORMED_ALT", record.filter)

class ExcludeMissingAltTestCase(test_case.JacquardBaseTestCase):
    def test_metaheader(self):
        self.assertEquals('##FILTER=<ID=JQ_EXCLUDE_MISSING_ALT,Description="The alternate allele is missing for this variant record.",Source="Jacquard",Version={}>'.format(__version__),
                          translate._ExcludeMissingAlt().metaheader)

    def test_add_tag_value_validAltNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "C", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_validIndelAltNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "ACGT", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_validIndelAltEdgecaseNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "ACGTNacgtn,*.", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_missingAltUpstreamDeletionNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "*", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_missingAltNullReplacesFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", ".", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("JQ_EXCLUDE_MISSING_ALT", record.filter)

    def test_add_tag_value_missingAltBothNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", ".*", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

class TranslateFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_translate(self):
        with TempDirectory() as output_file:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir,
                                          "functional_tests",
                                          "01_translate")
            input_file = os.path.join(module_testdir, "input")

            command = ["translate", input_file, output_file.path, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

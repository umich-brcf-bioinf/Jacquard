#pylint: disable=line-too-long,too-many-public-methods,invalid-name
#pylint: disable=missing-docstring,protected-access,too-few-public-methods
#pylint: disable=too-many-arguments,too-many-instance-attributes
from argparse import Namespace
import os

from testfixtures import TempDirectory

from jacquard import __version__, vcf
import jacquard.translate as translate
import test.test_case as test_case
from test.vcf_test import MockTranslatedVcfReader, MockVcfReader, MockCaller, MockTag, MockWriter

class TranslateTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(TranslateTestCase, self).setUp()

    def test_predict_output(self):
        with TempDirectory() as input_dir:
            input_dir.write("A.vcf", "##source=strelka\n#colHeader")
            input_dir.write("B.vcf", "##source=strelka\n#colHeader")
            args = Namespace(input=input_dir.path)

            desired_output_files = translate._predict_output(args)
            expected_desired_output_files = set(["A.translatedTags.vcf",
                                                 "B.translatedTags.vcf"])

            self.assertEquals(expected_desired_output_files, desired_output_files)

    def test_translate_files(self):
        mock_caller = MockCaller(metaheaders=["##mockCallerMetaheader1"])
        mock_tags = [MockTag(metaheader="##foo"), MockTag(metaheader="##bar")]
        reader = MockTranslatedVcfReader(MockVcfReader(), mock_tags, mock_caller)
        writer = MockWriter()
        execution_context = []
        translate._translate_files(reader, writer, execution_context)

        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertTrue(reader.add_tag_class_called)
        self.assertTrue(reader.vcf_records_called)

        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)

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
        record = vcf.VcfRecord("chr1", "42", "A", "ACGT*", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_validIndelAltEdgecaseNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "ACGTNacgtn,*.", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

    def test_add_tag_value_missingAltUpstreamDeletionReplacesFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", "*", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("JQ_EXCLUDE_MISSING_ALT", record.filter)

    def test_add_tag_value_missingAltNullReplacesFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", ".", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("JQ_EXCLUDE_MISSING_ALT", record.filter)

    def test_add_tag_value_missingAltBothNoFilter(self):
        record = vcf.VcfRecord("chr1", "42", "A", ".*", vcf_filter="PASS")
        translate._ExcludeMissingAlt().add_tag_values(record)
        self.assertEquals("PASS", record.filter)

class TranslateFunctionalTestCase(test_case.JacquardBaseTestCase):
    def Xtest_translate(self):
        with TempDirectory() as output_file:
            test_dir = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(test_dir,
                                          "functional_tests",
                                          "01_translate")
            input_file = os.path.join(module_testdir, "input")

            command = ["translate", input_file, output_file.path, "--force"]
            expected_dir = os.path.join(module_testdir, "benchmark")

            self.assertCommand(command, expected_dir)

#pylint:disable=line-too-long,too-many-public-methods
from __future__ import print_function, absolute_import, division

import jacquard.variant_caller_transforms.common_tags as common_tags
from jacquard.utils.vcf import VcfRecord
import test.utils.test_case as test_case
from jacquard.utils.utils import JQException


class JacquardTagTestCase(test_case.JacquardBaseTestCase):
    def test_allele_freq_tag(self):
        tag = common_tags.JacquardTag.ALLELE_FREQ_TAG
        self.assertEquals("AF", tag.abbreviation)
        self.assertEquals("Float", tag.vcf_type)
        self.assertEquals("A", tag.vcf_number)

    def test_depth_tag(self):
        tag = common_tags.JacquardTag.DEPTH_TAG
        self.assertEquals("DP", tag.abbreviation)
        self.assertEquals("Integer", tag.vcf_type)
        self.assertEquals("1", tag.vcf_number)

    def test_genotype_tag(self):
        tag = common_tags.JacquardTag.GENOTYPE_TAG
        self.assertEquals("GT", tag.abbreviation)
        self.assertEquals("String", tag.vcf_type)
        self.assertEquals("1", tag.vcf_number)

    def test_somatic_tag(self):
        tag = common_tags.JacquardTag.SOMATIC_TAG
        self.assertEquals("HC_SOM", tag.abbreviation)
        self.assertEquals("Integer", tag.vcf_type)
        self.assertEquals("1", tag.vcf_number)

    def test_add_tag_values_raisesNotImplementedError(self):
        class FakeTag(common_tags.JacquardTag):
            def __init__(self): pass
        tag = FakeTag()
        self.assertRaises(NotImplementedError,
                          tag.add_tag_values,
                          VcfRecord("1", "42", "A", "C")
                          )

    def test_metaheader(self):
        tag_type = common_tags.JacquardTag.GENOTYPE_TAG
        class MyTag(common_tags.JacquardTag):
            def add_tag_values(self, vcf_record): pass
        actual_tag = MyTag("SK", tag_type, "foo bar baz")
        expected_metaheader = ('##FORMAT=<ID=JQ_SK_GT,Number=1,Type=String,'
                               'Description="foo bar baz">')
        self.assertEquals(expected_metaheader, actual_tag.metaheader)

    def test_metaheader_raisesExceptionIfEmbeddedQuotesInDescription(self):
        tag_type = common_tags.JacquardTag.GENOTYPE_TAG
        class MyTag(common_tags.JacquardTag):
            def add_tag_values(self, vcf_record): pass
        MyTag("SK", tag_type, "Single quotes are 'ok'")
        self.assertRaisesRegexp(JQException,
                                r'Metaheader descriptions cannot contain double quotes: \[Double quotes are "not ok"\]',
                                MyTag,
                                "SK",
                                tag_type,
                                'Double quotes are "not ok"')

    def test_tag_id(self):
        tag_type = common_tags.JacquardTag.GENOTYPE_TAG
        class MyTag(common_tags.JacquardTag):
            def add_tag_values(self, vcf_record): pass
        actual_tag = MyTag("SK", tag_type, "foo bar baz")
        self.assertEquals("JQ_SK_GT", actual_tag.tag_id)


class ReportedTagTestCase(test_case.JacquardBaseTestCase):
    def test_reported_tag_metaheader(self):
        reported_tag = common_tags.ReportedTag("foo_")
        self.assertEquals(('##FORMAT=<ID={}{},'
                           'Number=1,'
                           'Type=Integer,'
                           'Description="1 = variant present in original VCF">')\
                          .format("foo_",
                                  common_tags.CALLER_REPORTED_TAG),
                          reported_tag.metaheader)

    def test_reported_tag_format(self):
        reported_tag = common_tags.ReportedTag("foo_")
        actual_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2|SA.1:SA.2|SB.1:SB.2\n")
        actual_vcf_record = VcfRecord.parse_record(actual_line, ["SA", "SB"])
        expected_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:foo_" +
                                   common_tags.CALLER_REPORTED_TAG +
                                   "|SA.1:SA.2:1|SB.1:SB.2\n")
        expected_vcf_record = VcfRecord.parse_record(expected_line, ["SA", "SB"])
        reported_tag.add_tag_values(actual_vcf_record)
        self.assertEquals(expected_vcf_record, actual_vcf_record)

    def test_passed_tag_metaheader(self):
        passed_tag = common_tags.PassedTag("foo_")
        self.assertEquals(('##FORMAT=<ID={}{},'
                           'Number=1,Type=Integer,'
                           'Description="1 = variant FILTER is PASS in '
                           'original VCF">').format("foo_",
                                                    common_tags.CALLER_PASSED_TAG),
                          passed_tag.metaheader)

    def test_passed_tag_format(self):
        passed_tag = common_tags.PassedTag("foo_")
        actual_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2|SA.1:SA.2|SB.1:SB.2\n")
        actual_vcf_record = VcfRecord.parse_record(actual_line, ["SA", "SB"])
        expected_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:foo_" +
                                   common_tags.CALLER_PASSED_TAG +
                                   "|SA.1:SA.2:1|SB.1:SB.2\n")
        expected_vcf_record = VcfRecord.parse_record(expected_line, ["SA", "SB"])
        passed_tag.add_tag_values(actual_vcf_record)
        self.assertEquals(expected_vcf_record, actual_vcf_record)

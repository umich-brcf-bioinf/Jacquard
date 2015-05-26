#pylint:disable=line-too-long,too-many-public-methods,too-few-public-methods
#pylint:disable=invalid-name, abstract-method, multiple-statements,
#pylint:disable=super-init-not-called
from __future__ import print_function, absolute_import, division

from collections import OrderedDict

from jacquard.utils.utils import JQException
from jacquard.utils.vcf import VcfRecord
import jacquard.variant_caller_transforms.common_tags as common_tags
import test.utils.test_case as test_case


class AbstractJacquardTagTestCase(test_case.JacquardBaseTestCase):
    def test_allele_freq_tag(self):
        tag = common_tags.ALLELE_FREQ_TAG
        self.assertEquals("AF", tag.abbreviation)
        self.assertEquals("Float", tag.vcf_type)
        self.assertEquals("A", tag.vcf_number)

    def test_depth_tag(self):
        tag = common_tags.DEPTH_TAG
        self.assertEquals("DP", tag.abbreviation)
        self.assertEquals("Integer", tag.vcf_type)
        self.assertEquals("1", tag.vcf_number)

    def test_genotype_tag(self):
        tag = common_tags.GENOTYPE_TAG
        self.assertEquals("GT", tag.abbreviation)
        self.assertEquals("String", tag.vcf_type)
        self.assertEquals("1", tag.vcf_number)

    def test_somatic_tag(self):
        tag = common_tags.SOMATIC_TAG
        self.assertEquals("HC_SOM", tag.abbreviation)
        self.assertEquals("Integer", tag.vcf_type)
        self.assertEquals("1", tag.vcf_number)

    def test_add_tag_values_raisesNotImplementedError(self):
        class FakeTag(common_tags.AbstractJacquardTag):
            def __init__(self): pass
        tag = FakeTag()
        self.assertRaises(NotImplementedError,
                          tag.add_tag_values,
                          VcfRecord("1", "42", "A", "C")
                          )

    def test_metaheader(self):
        tag_type = common_tags.GENOTYPE_TAG
        class MyTag(common_tags.AbstractJacquardTag):
            def add_tag_values(self, vcf_record): pass
        actual_tag = MyTag("SK", tag_type, "foo bar baz")
        expected_metaheader = ('##FORMAT=<ID=JQ_SK_GT,Number=1,Type=String,'
                               'Description="foo bar baz">')
        self.assertEquals(expected_metaheader, actual_tag.metaheader)

    def test_metaheader_raisesExceptionIfEmbeddedQuotesInDescription(self):
        tag_type = common_tags.GENOTYPE_TAG
        class MyTag(common_tags.AbstractJacquardTag):
            def add_tag_values(self, vcf_record): pass
        MyTag("SK", tag_type, "Single quotes are 'ok'")
        self.assertRaisesRegexp(JQException,
                                r'Metaheader descriptions cannot contain double quotes: \[Double quotes are "not ok"\]',
                                MyTag,
                                "SK",
                                tag_type,
                                'Double quotes are "not ok"')

    def test_tag_id(self):
        tag_type = common_tags.GENOTYPE_TAG
        class MyTag(common_tags.AbstractJacquardTag):
            def add_tag_values(self, vcf_record): pass
        actual_tag = MyTag("SK", tag_type, "foo bar baz")
        self.assertEquals("JQ_SK_GT", actual_tag.tag_id)

    def test_get_matching_tags_exact(self):
        format_tags = OrderedDict(sorted({"JQ_DP": "32", "JQ_GT": "0/1", "JQ_AF": "0.2"}.items()))
        tag_type = common_tags.TagType("GT", "String", "1")
        actual_tags = common_tags.AbstractJacquardTag.get_matching_tags(format_tags,
                                                                        tag_type)

        expected_tags = OrderedDict({"JQ_GT": "0/1"})

        self.assertEquals(expected_tags, actual_tags)

    def test_get_matching_tags_regex(self):
        format_tags = OrderedDict(sorted({"JQ_MT_DP": "32", "JQ_VS_DP": "43", "JQ_MT_GT": "0/1", "JQ_VS_GT": "0/0"}.items()))
        tag_type = common_tags.TagType("GT", "String", "1")
        actual_tags = common_tags.AbstractJacquardTag.get_matching_tags(format_tags,
                                                                        tag_type)

        expected_tags = OrderedDict(sorted({"JQ_MT_GT": "0/1", "JQ_VS_GT": "0/0"}.items()))
        self.assertEquals(expected_tags, actual_tags)
        self.assertEquals(["JQ_MT_GT", "JQ_VS_GT"], actual_tags.keys())

    def test_get_matching_tags_regexSuffix(self):
        format_tags = OrderedDict(sorted({"JQ_MT_DP": "32", "JQ_GT_FOO": "bar", "JQ_MT_GT": "0/1", "JQ_VS_GT": "0/0"}.items()))
        tag_type = common_tags.TagType("GT", "String", "1")
        actual_tags = common_tags.AbstractJacquardTag.get_matching_tags(format_tags,
                                                                        tag_type)

        expected_tags = OrderedDict(sorted({"JQ_MT_GT": "0/1", "JQ_VS_GT": "0/0"}.items()))
        self.assertEquals(expected_tags, actual_tags)
        self.assertEquals(["JQ_MT_GT", "JQ_VS_GT"], actual_tags.keys())

    def test_get_matching_tags_regexOnlyJacquard(self):
        format_tags = OrderedDict(sorted({"JQ_SK_GT": "0/1", "GT": "0/1"}.items()))
        tag_type = common_tags.TagType("GT", "String", "1")
        actual_tags = common_tags.AbstractJacquardTag.get_matching_tags(format_tags,
                                                                        tag_type)

        expected_tags = OrderedDict({"JQ_SK_GT": "0/1"})
        self.assertEquals(expected_tags, actual_tags)

    def test_get_matching_caller_abbrevs(self):
        format_tags = OrderedDict(sorted({"JQ_VS_DP": "32", "JQ_VS_GT": "0/1", "JQ_VS_AF": "0.2"}.items()))
        tag_type = common_tags.TagType("GT", "String", "1")
        actual_tags = common_tags.AbstractJacquardTag.get_matching_caller_abbrevs(format_tags,
                                                                                  tag_type)

        expected_tags = OrderedDict({"VS": "0/1"})
        self.assertEquals(expected_tags, actual_tags)

    def test_get_matching_caller_abbrevs_regex(self):
        format_tags = OrderedDict(sorted({"JQ_MT_DP": "32", "JQ_VS_DP": "43", "JQ_MT_GT": "0/1", "JQ_VS_GT": "0/0"}.items()))
        tag_type = common_tags.TagType("GT", "String", "1")
        actual_tags = common_tags.AbstractJacquardTag.get_matching_caller_abbrevs(format_tags,
                                                                                  tag_type)

        expected_tags = OrderedDict(sorted({"MT": "0/1", "VS": "0/0"}.items()))
        self.assertEquals(expected_tags, actual_tags)

    def test_get_matching_caller_abbrevs_regexSuffix(self):
        format_tags = OrderedDict(sorted({"JQ_MT_DP": "32", "JQ_GT_FOO": "bar", "JQ_MT_GT": "0/1", "JQ_VS_GT": "0/0"}.items()))
        tag_type = common_tags.TagType("GT", "String", "1")
        actual_tags = common_tags.AbstractJacquardTag.get_matching_caller_abbrevs(format_tags,
                                                                                  tag_type)

        expected_tags = OrderedDict(sorted({"MT": "0/1", "VS": "0/0"}.items()))
        self.assertEquals(expected_tags, actual_tags)

    def test_get_matching_caller_abbrevs_regexOnlyJacquard(self):
        format_tags = OrderedDict(sorted({"JQ_SK_GT": "0/1", "GT": "0/1"}.items()))
        tag_type = common_tags.TagType("GT", "String", "1")
        actual_tags = common_tags.AbstractJacquardTag.get_matching_caller_abbrevs(format_tags,
                                                                                  tag_type)

        expected_tags = OrderedDict({"SK": "0/1"})
        self.assertEquals(expected_tags, actual_tags)

class ReportedTagTestCase(test_case.JacquardBaseTestCase):
    def test_reported_tag_metaheader(self):
        reported_tag = common_tags.ReportedTag("foo")
        self.assertEquals(('##FORMAT=<ID={},'
                           'Number=1,'
                           'Type=Integer,'
                           'Description="1 = variant present in original VCF">')\
                          .format("JQ_foo_CALLER_REPORTED"),
                          reported_tag.metaheader)

    def test_reported_tag_format(self):
        reported_tag = common_tags.ReportedTag("foo_")
        actual_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2|SA.1:SA.2|SB.1:SB.2\n")
        actual_vcf_record = VcfRecord.parse_record(actual_line, ["SA", "SB"])
        expected_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:foo_CALLER_REPORTED"
                                   "|SA.1:SA.2:1|SB.1:SB.2\n")
        expected_vcf_record = VcfRecord.parse_record(expected_line, ["SA", "SB"])
        reported_tag.add_tag_values(actual_vcf_record)
        self.assertEquals(expected_vcf_record, actual_vcf_record)

    def test_passed_tag_metaheader(self):
        passed_tag = common_tags.PassedTag("foo")
        self.assertEquals(('##FORMAT=<ID={},'
                           'Number=1,Type=Integer,'
                           'Description="1 = variant FILTER is PASS in '
                           'original VCF">').format("JQ_foo_CALLER_PASSED"),
                          passed_tag.metaheader)

    def test_passed_tag_format(self):
        passed_tag = common_tags.PassedTag("foo")
        actual_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2|SA.1:SA.2|SB.1:SB.2\n")
        actual_vcf_record = VcfRecord.parse_record(actual_line, ["SA", "SB"])
        expected_line = self.entab("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:foo_CALLER_REPORTED"
                                   "|SA.1:SA.2:1|SB.1:SB.2\n")
        expected_vcf_record = VcfRecord.parse_record(expected_line, ["SA", "SB"])
        passed_tag.add_tag_values(actual_vcf_record)
        self.assertEquals(expected_vcf_record, actual_vcf_record)

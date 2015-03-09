#pylint:disable=line-too-long,too-many-public-methods
from jacquard import __version__
from jacquard.vcf import VcfRecord
import jacquard.variant_callers.common_tags as common_tags
import test.test_case as test_case

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

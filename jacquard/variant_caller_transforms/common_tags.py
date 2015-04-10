"""Common tags used by several callers."""
from __future__ import print_function, absolute_import, division


CALLER_REPORTED_TAG = "CALLER_REPORTED"
CALLER_PASSED_TAG = "CALLER_PASSED"

class ReportedTag(object):
    """Tracks whether the caller reported this variant (i.e. it's in the VCF).

    This tag could be inferred through the presence of other tags, but adding
    it explicitly simplifies how summary tags are generated.
    """
    #pylint: disable=too-few-public-methods
    def __init__(self, tag_name):
        self.tag_name = tag_name
        self.metaheader = ('##FORMAT=<ID={}{},'
                           'Number=1,'
                           'Type=Integer,'
                           #pylint: disable=line-too-long
                           'Description="1 = variant present in original VCF">')\
                           .format(self.tag_name,
                                   CALLER_REPORTED_TAG)

    def add_tag_values(self, vcf_record):
        sample_values = {}
        for sample in vcf_record.sample_tag_values:
            sample_values[sample] = 1
        vcf_record.add_sample_tag_value(self.tag_name+CALLER_REPORTED_TAG,
                                        sample_values)

class PassedTag(object):
    """Tracks whether the caller passed this variant.

    This enables a useful summary tag. Note that Jaquard may flag a variant
    that originally passed as invalid (e.g. invalid Varscan ALT value). In
    this case, this tag represents what the origin caller thought and not
    Jacquard's opinion.
    """
    #pylint: disable=too-few-public-methods
    def __init__(self, tag_name):
        self.tag_name = tag_name
        self.metaheader = ('##FORMAT=<ID={}{},'
                           'Number=1,Type=Integer,'
                           'Description="1 = variant FILTER is PASS in '
                           'original VCF">').format(self.tag_name,
                                                    CALLER_PASSED_TAG)

    def add_tag_values(self, vcf_record):
        sample_values = {}
        for sample in vcf_record.sample_tag_values:
            if vcf_record.filter == "PASS":
                sample_values[sample] = 1
            else:
                sample_values[sample] = 0
        vcf_record.add_sample_tag_value(self.tag_name + CALLER_PASSED_TAG,
                                        sample_values)

"""Common tags used by several callers."""
from __future__ import print_function, absolute_import, division
import abc
import jacquard.utils.utils as utils

CALLER_REPORTED_TAG = "CALLER_REPORTED"
CALLER_PASSED_TAG = "CALLER_PASSED"

class AbstractJacquardTag(object):
    #pylint: disable=abstract-class-not-used, too-few-public-methods
    __metaclass__ = abc.ABCMeta

    FORMAT = ('##FORMAT=<ID={},Number={},'
              'Type={},Description="{}">')

    class _TagType(object):
        #pylint: disable=too-few-public-methods
        def __init__(self, abbreviation, vcf_type, vcf_number):
            self.abbreviation = abbreviation
            self.vcf_type = vcf_type
            self.vcf_number = vcf_number

    DEPTH_TAG = _TagType("DP", "Integer", "1")
    GENOTYPE_TAG = _TagType("GT", "String", "1")
    ALLELE_FREQ_TAG = _TagType("AF", "Float", "A")
    SOMATIC_TAG = _TagType("HC_SOM", "Integer", "1")
    CALLER_REPORTED_TAG = _TagType("CALLER_REPORTED", "Integer", "1")
    CALLER_PASSED_TAG = _TagType("CALLER_PASSED", "Integer", "1")

    def __init__(self, variant_caller_abbrev, tag_type, description):
        if '"' in description:
            raise utils.JQException(("Metaheader descriptions cannot contain "
                                    "double quotes: [{}]"),
                                    description)
        self.tag_id = "JQ_{}_{}".format(variant_caller_abbrev,
                                        tag_type.abbreviation)
        self.metaheader = AbstractJacquardTag.FORMAT.format(self.tag_id,
                                                            tag_type.vcf_number,
                                                            tag_type.vcf_type,
                                                            description)

    def add_tag_values(self, vcf_record):
        raise NotImplementedError()


class ReportedTag(AbstractJacquardTag):
    """Tracks whether the caller reported this variant (i.e. it's in the VCF).

    This tag could be inferred through the presence of other tags, but adding
    it explicitly simplifies how summary tags are generated.
    """
    #pylint: disable=too-few-public-methods
    
    def __init__(self, caller_abbreviation):
        super(self.__class__,
              self).__init__(caller_abbreviation,
                             AbstractJacquardTag.CALLER_REPORTED_TAG,
                             '1 = variant present in original VCF')

    def add_tag_values(self, vcf_record):
        sample_values = {}
        for sample in vcf_record.sample_tag_values:
            sample_values[sample] = 1
        vcf_record.add_sample_tag_value(self.tag_id,
                                        sample_values)

class PassedTag(AbstractJacquardTag):
    """Tracks whether the caller passed this variant.

    This enables a useful summary tag. Note that Jaquard may flag a variant
    that originally passed as invalid (e.g. invalid Varscan ALT value). In
    this case, this tag represents what the origin caller thought and not
    Jacquard's opinion.
    """
    #pylint: disable=too-few-public-methods
    def __init__(self, caller_abbreviation):
        super(self.__class__,
              self).__init__(caller_abbreviation,
                             AbstractJacquardTag.CALLER_PASSED_TAG,
                             '1 = variant FILTER is PASS in original VCF')

    def add_tag_values(self, vcf_record):
        sample_values = {}
        for sample in vcf_record.sample_tag_values:
            if vcf_record.filter == "PASS":
                sample_values[sample] = 1
            else:
                sample_values[sample] = 0
        vcf_record.add_sample_tag_value(self.tag_id,
                                        sample_values)

"""Common tags used by several callers."""
from __future__ import print_function, absolute_import, division

import abc
import re
from collections import OrderedDict

import jacquard.utils.utils as utils


CALLER_REPORTED = "CALLER_REPORTED"
CALLER_PASSED = "CALLER_PASSED"

class TagType(object):
    #pylint: disable=too-few-public-methods
    def __init__(self,
                 abbreviation,
                 vcf_type,
                 vcf_number,
                 metaheader_type='FORMAT'):

        self.abbreviation = abbreviation
        self.vcf_type = vcf_type
        self.vcf_number = vcf_number
        self.metaheader_type = metaheader_type

DEPTH_TAG = TagType("DP", "Integer", "1")
GENOTYPE_TAG = TagType("GT", "String", "1")
ALLELE_FREQ_TAG = TagType("AF", "Float", "A")
SOMATIC_TAG = TagType("HC_SOM", "Integer", "1")

CALLER_REPORTED_TAG = TagType("CALLER_REPORTED", "Integer", "1")
CALLER_PASSED_TAG = TagType("CALLER_PASSED", "Integer", "1")

class AbstractJacquardTag(object):
    #pylint: disable=abstract-class-not-used, too-few-public-methods
    __metaclass__ = abc.ABCMeta

    FORMAT = ('##{}=<ID={},Number={},'
              'Type={},Description="{}">')

    def __init__(self, variant_caller_abbrev, tag_type, description):

        if '"' in description:
            raise utils.JQException(("Metaheader descriptions cannot contain "
                                    "double quotes: [{}]"),
                                    description)
        self.tag_type = tag_type
        self.tag_id = "JQ_{}_{}".format(variant_caller_abbrev,
                                        tag_type.abbreviation)
        self.metaheader = AbstractJacquardTag.FORMAT.format(self.tag_type\
                                                            .metaheader_type,
                                                            self.tag_id,
                                                            self.tag_type\
                                                            .vcf_number,
                                                            self.tag_type\
                                                            .vcf_type,
                                                            description)

    def add_tag_values(self, vcf_record):
        raise NotImplementedError()

    @staticmethod
    def get_matching_tags(format_tags, tag_type):
        matching_format_tags = OrderedDict()
        for key, value in list(format_tags.items()):
            regex = re.compile(r"^JQ.*{}$".format(tag_type.abbreviation))
            if re.search(regex, key):
                matching_format_tags[key] = value

        return matching_format_tags

    @staticmethod
    def get_matching_caller_abbrevs(format_tags, tag_type):
        matching_format_tags = OrderedDict()
        for key, value in list(format_tags.items()):
            regex = re.compile(r"^JQ_(.*)_{}$".format(tag_type.abbreviation))
            match = re.search(regex, key)
            if match:
                matching_format_tags[match.group(1)] = value

        return matching_format_tags

class ReportedTag(AbstractJacquardTag):
    """Tracks whether the caller reported this variant (i.e. it's in the VCF).

    This tag could be inferred through the presence of other tags, but adding
    it explicitly simplifies how summary tags are generated.
    """
    #pylint: disable=too-few-public-methods

    def __init__(self, caller_abbreviation):
        super(self.__class__,
              self).__init__(caller_abbreviation,
                             CALLER_REPORTED_TAG,
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
                             CALLER_PASSED_TAG,
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

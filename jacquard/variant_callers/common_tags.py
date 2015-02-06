#pylint: disable=too-few-public-methods, unused-argument
from __future__ import print_function, absolute_import
import jacquard.utils as utils

CALLER_REPORTED_TAG = "CALLER_REPORTED"
CALLER_PASSED_TAG = "CALLER_PASSED"

class ReportedTag(object):
    def __init__(self, tag_name):
        self.tag_name = tag_name
        self.metaheader = ('##FORMAT=<ID={}{},'
                           'Number=1,Type=Integer,'
                           'Description="0 = variant absent from original VCF; '
                           '1 = variant present in original VCF",'
                           'Source="Jacquard",'
                           'Version="{}">').format(self.tag_name,
                                                   CALLER_REPORTED_TAG,
                                                   utils.__version__)

    def add_tag_values(self, vcf_record):
        sample_values = {}
        for sample in vcf_record.sample_tag_values:
            sample_values[sample] = 1
        vcf_record.add_sample_tag_value(self.tag_name+CALLER_REPORTED_TAG,
                                        sample_values)

class PassedTag(object):
    def __init__(self, tag_name):
        self.tag_name = tag_name
        self.metaheader = ('##FORMAT=<ID={}{},'
                           'Number=1,Type=Integer,'
                           'Description="0 = variant FILTER is not PASS in '
                           'original VCF; '
                           '1 = variant FILTER is PASS in original VCF",'
                           'Version="{}">').format(self.tag_name,
                                                   CALLER_PASSED_TAG,
                                                   utils.__version__)

    def add_tag_values(self, vcf_record):
        sample_values = {}
        for sample in vcf_record.sample_tag_values:
            if vcf_record.filter == "PASS":
                sample_values[sample] = 1
            else:
                sample_values[sample] = 0
        vcf_record.add_sample_tag_value(self.tag_name + CALLER_PASSED_TAG,
                                        sample_values)

"""Classes to summarize data for a sample-variant or variant as a whole.

A collection of individual Tag classes hold the metaheader and logic to
transform Jacquard-standardized VcfRecords.
"""
from __future__ import print_function, absolute_import, division

from collections import defaultdict
from decimal import Decimal

import jacquard.utils.utils as utils
import jacquard.variant_caller_transforms.common_tags as common_tags
from  jacquard.variant_caller_transforms.common_tags import TagType as TagType

SUMMARY_TAG = "SUMMARY"

JQ_SUMMARY_TAG = "JQ_SUMMARY_"
JQ_REPORTED = "CALLERS_REPORTED_COUNT"
JQ_SAMPLES_REPORTED = "SAMPLES_REPORTED_COUNT"
JQ_PASSED = "CALLERS_PASSED_COUNT"
JQ_SAMPLES_PASSED = "SAMPLES_PASSED_COUNT"

CALLERS_REPORTED_LIST = TagType("CALLERS_REPORTED_LIST", "String", ".")
CALLERS_PASSED_LIST = TagType("CALLERS_PASSED_LIST", "String", ".")
CALLERS_REPORTED_COUNT = TagType("CALLERS_REPORTED_COUNT", "Integer", "1")
CALLERS_PASSED_COUNT = TagType("CALLERS_PASSED_COUNT", "Integer", "1")

SAMPLES_REPORTED = TagType("SAMPLES_REPORTED_COUNT", "Integer", "1", "INFO")
SAMPLES_PASSED = TagType("SAMPLES_PASSED_COUNT", "Integer", "1", "INFO")

SUMMARY_GENOTYPE = TagType("HC_GT", "String", "1")
SUMMARY_ALLELE_FREQ_RANGE = TagType("AF_RANGE", "Float", "1")
SUMMARY_ALLELE_FREQ_AVG = TagType("AF_AVERAGE", "Float", "1")
SUMMARY_DEPTH_RANGE = TagType("DP_RANGE", "Float", "1")
SUMMARY_DEPTH_AVG = TagType("DP_AVERAGE", "Float", "1")
SUMMARY_SOMATIC_COUNT = TagType("SOM_COUNT", "Integer", "1")

def _aggregate_numeric_values(values, function):
    split_values = [x.split(",") for x in values]
    lengths = set([len(x) for x in split_values])
    if len(lengths) > 1:
        return 0

    #pylint: disable=star-args
    transposed_values = [list(x) for x in zip(*split_values)]
    string_values = []
    for group in transposed_values:
        string_values.append(function([Decimal(i) for i in group]))

    return ",".join(string_values)

def _get_non_null_values(record, sample, tag_type):
    values = []
    try:
        desired_tags = common_tags.AbstractJacquardTag\
                       .get_matching_tags(record.sample_tag_values[sample],
                                          tag_type)
        for value in list(desired_tags.values()):
            if value != ".":
                values.append(value)
    except KeyError:
        raise utils.JQException("Sample [{}] was not recognized".format(sample))

    return values

def _average(numeric_values):
    average = str(sum(numeric_values)/Decimal(len(numeric_values)))
    return utils.round_digits(average)

def _range(numeric_values):
    if len(numeric_values) > 1:
        value_range = str(max(numeric_values) - min(numeric_values))
        return utils.round_digits(value_range)
    return "."

def _count(numeric_values):
    return str(sum(i == 1 for i in numeric_values))

def _build_new_tags(vcf_record, tags, sample):
    desired_tags = []
    for tag in tags:
        desired_tag = vcf_record.sample_tag_values[sample][tag].split(",")

        #don't include null values in avg calculation
        altered_tag = [x for x in desired_tag if x != "."]
        if len(altered_tag) != 0:
            desired_tags.append(altered_tag)

    return desired_tags

def _add_caller_list_values(vcf_record, tag_type, tag_id):
    sample_tag = {}
    for sample, tags in vcf_record.sample_tag_values.items():
        sample_tag[sample] = []
        desired_tags = common_tags.AbstractJacquardTag\
                       .get_matching_caller_abbrevs(tags, tag_type)
        for tag in desired_tags:
            if desired_tags[tag] != "0" and desired_tags[tag] != ".":
                sample_tag[sample].append(tag)

    for sample in sample_tag:
        sample_tag[sample] = ",".join(sorted(sample_tag[sample]))
        if not sample_tag[sample]:
            sample_tag[sample] = "."
    vcf_record.add_sample_tag_value(tag_id, sample_tag)

def _add_caller_count_values(vcf_record, tag_type, tag_id):
    sample_tag = {}
    for sample, all_tags in vcf_record.sample_tag_values.items():
        sample_tag[sample] = 0
        tags = common_tags.AbstractJacquardTag.get_matching_tags(all_tags,
                                                                 tag_type)
        for tag in tags:
            if tags[tag] != "0" and tags[tag] != ".":
                sample_tag[sample] += int(tags[tag])
    vcf_record.add_sample_tag_value(tag_id, sample_tag)

def _add_sample_count_values(vcf_record,
                             jq_callers_global_variable,
                             tag_id):
    count = 0
    for tags in vcf_record.sample_tag_values.values():
        tag_key = JQ_SUMMARY_TAG + jq_callers_global_variable
        if tag_key in tags and tags[tag_key] != "0" and tags[tag_key] != ".":
            count += 1

    vcf_record.add_info_field("=".join([tag_id, str(count)]))

class _CallersReportedListTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods

    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             CALLERS_REPORTED_LIST,
                             ('Comma-separated list variant callers '
                             'which listed this variant in the Jacquard tagged '
                             'VCF'))

    def add_tag_values(self, vcf_record):
        _add_caller_list_values(vcf_record,
                                common_tags.CALLER_REPORTED_TAG,
                                self.tag_id)

class _CallersReportedTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             CALLERS_REPORTED_COUNT,
                             ('Count of variant callers which listed this '
                              'variant in the Jacquard tagged VCF'))

    def add_tag_values(self, vcf_record):
        _add_caller_count_values(vcf_record,
                                 common_tags.CALLER_REPORTED_TAG,
                                 self.tag_id)

class _CallersPassedListTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods

    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             CALLERS_PASSED_LIST,
                             ('Comma-separated list of variant caller short-'
                              'names where FILTER = PASS for this variant in '
                              'the Jacquard tagged VCF'))

    def add_tag_values(self, vcf_record):
        _add_caller_list_values(vcf_record,
                                common_tags.CALLER_PASSED_TAG,
                                self.tag_id)

class _CallersPassedTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods

    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             CALLERS_PASSED_COUNT,
                             ('Count of variant callers where FILTER = PASS'
                              ' for this variant in the Jacquard tagged VCF'))

    def add_tag_values(self, vcf_record):
        _add_caller_count_values(vcf_record,
                                 common_tags.CALLER_PASSED_TAG,
                                 self.tag_id)

class _SamplesReported(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             SAMPLES_REPORTED,
                             ('Count of samples where this variant appeared in '
                              'any of the Jacquard tagged VCFs (regardless of '
                              'quality/filtering)'))

    def add_tag_values(self, vcf_record):
        _add_sample_count_values(vcf_record, JQ_REPORTED, self.tag_id)

class _SamplesPassed(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             SAMPLES_PASSED,
                             ('Count of samples where a variant caller passed '
                              'the filter in any of the Jacquard tagged VCFs'))

    def add_tag_values(self, vcf_record):
        _add_sample_count_values(vcf_record, JQ_PASSED, self.tag_id)

class _HCGenotypeTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods

    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             SUMMARY_GENOTYPE,
                             ('High confidence consensus genotype (inferred '
                              'from JQ_*_GT and JQ_*_CALLER_PASSED). Majority '
                              'rules; ties go to the least unusual variant '
                              '(0/1>0/2>1/1). Variants which failed their '
                              'filter are ignored. Phasing is removed.'))

    @staticmethod
    def _prioritize_genotype(gts):

        def _dephase(gt):
            result = gt
            if '|' in gt:
                result = '/'.join(sorted(gt.split('|')))
            return result

        def _break_ties(value):
            if value == "0/0":
                return "99"
            return value

        dephased_gts = list(map(_dephase, gts))
        count = defaultdict(int)
        for gt in dephased_gts:
            count[gt] += 1

        sorted_gts = sorted(dephased_gts,
                            key=lambda x: (-count[x], _break_ties(x)))
        return sorted_gts[0]

    def add_tag_values(self, vcf_record):
        new_sample_tag_values = {}
        for sample in vcf_record.sample_tag_values:
            tag_values = _get_non_null_values(vcf_record,
                                              sample,
                                              common_tags.GENOTYPE_TAG)
            genotype = "."
            if tag_values:
                genotype = self._prioritize_genotype(tag_values)
            new_sample_tag_values[sample] = genotype

        vcf_record.add_sample_tag_value(self.tag_id,
                                        new_sample_tag_values)

class _AlleleFreqRangeTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods

    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             SUMMARY_ALLELE_FREQ_RANGE,
                             ('Max(allele frequency) - min (allele frequency) '
                              'across recognized callers.'))

    def add_tag_values(self, record):
        new_sample_tag_values = {}
        for sample in record.sample_tag_values:
            tag_values = _get_non_null_values(record,
                                              sample,
                                              common_tags.ALLELE_FREQ_TAG)

            aggregated_values = "."
            if tag_values:
                aggregated_values = _aggregate_numeric_values(tag_values,
                                                              _range)
            if not aggregated_values:
                msg = "Error summarizing values {} at record [{}:{} {}]"
                raise utils.JQException(msg,
                                        list(tag_values),
                                        record.chrom,
                                        record.pos,
                                        sample)

            new_sample_tag_values[sample] = aggregated_values

        record.add_sample_tag_value(self.tag_id,
                                    new_sample_tag_values)

class _AlleleFreqAverageTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods

    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             SUMMARY_ALLELE_FREQ_AVG,
                             ('Average allele frequency across recognized '
                              'variant callers that reported frequency for '
                              'this sample-locus [average(JQ_*_AF)].'))

    def add_tag_values(self, record):
        new_sample_tag_values = {}
        for sample in record.sample_tag_values:
            tag_values = _get_non_null_values(record,
                                              sample,
                                              common_tags.ALLELE_FREQ_TAG)

            aggregated_values = "."
            if tag_values:
                aggregated_values = _aggregate_numeric_values(tag_values,
                                                              _average)
            if not aggregated_values:
                msg = "Error summarizing values {} at record [{}:{} {}]"
                raise utils.JQException(msg,
                                        list(tag_values),
                                        record.chrom,
                                        record.pos,
                                        sample)

            new_sample_tag_values[sample] = aggregated_values

        record.add_sample_tag_value(self.tag_id,
                                    new_sample_tag_values)

class _DepthRangeTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods

    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             SUMMARY_DEPTH_RANGE,
                             ('Max(depth) - min (depth) across recognized '
                              'callers.'))

    def add_tag_values(self, record):
        new_sample_tag_values = {}
        for sample in record.sample_tag_values:
            tag_values = _get_non_null_values(record,
                                              sample,
                                              common_tags.DEPTH_TAG)

            aggregated_values = "."
            if tag_values:
                aggregated_values = _aggregate_numeric_values(tag_values,
                                                              _range)
            if not aggregated_values:
                msg = "Error summarizing values {} at record [{}:{} {}]"
                raise utils.JQException(msg,
                                        list(tag_values),
                                        record.chrom,
                                        record.pos,
                                        sample)

            new_sample_tag_values[sample] = aggregated_values

        record.add_sample_tag_value(self.tag_id,
                                    new_sample_tag_values)

class _DepthAverageTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods

    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             SUMMARY_DEPTH_AVG,
                             ('Average depth across recognized '
                              'variant callers that reported depth for '
                              'this sample-locus; rounded to integer '
                              '[round(average(JQ_*_DP))].'))

    def add_tag_values(self, record):
        new_sample_tag_values = {}
        for sample in record.sample_tag_values:
            tag_values = _get_non_null_values(record,
                                              sample,
                                              common_tags.DEPTH_TAG)

            aggregated_values = "."
            if tag_values:
                aggregated_values = _aggregate_numeric_values(tag_values,
                                                              _average)
            if not aggregated_values:
                msg = "Error summarizing values {} at record [{}:{} {}]"
                raise utils.JQException(msg,
                                        list(tag_values),
                                        record.chrom,
                                        record.pos,
                                        sample)

            new_sample_tag_values[sample] = aggregated_values

        record.add_sample_tag_value(self.tag_id,
                                    new_sample_tag_values)

class _SomaticTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods

    def __init__(self):
        super(self.__class__,
              self).__init__(SUMMARY_TAG,
                             SUMMARY_SOMATIC_COUNT,
                             ('Count of recognized variant callers that '
                              'reported confident somatic call for this '
                              'sample-locus.'))

    def add_tag_values(self, record):
        new_sample_tag_values = {}
        for sample in record.sample_tag_values:
            tag_values = _get_non_null_values(record,
                                              sample,
                                              common_tags.SOMATIC_TAG)
            aggregated_values = "."
            if tag_values:
                aggregated_values = _aggregate_numeric_values(tag_values,
                                                              _count)
            if not aggregated_values:
                msg ="Error summarizing values {} at record [{}:{} {}]"
                raise utils.JQException(msg,
                                        list(tag_values),
                                        record.chrom,
                                        record.pos,
                                        sample)

            new_sample_tag_values[sample] = aggregated_values

        record.add_sample_tag_value(self.tag_id,
                                    new_sample_tag_values)

class SummarizeCaller(object):
    """Provides metaheaders for VcfReader; adds summary tags to VcfRecords."""
    def __init__(self):
        self.tags = [_CallersReportedTag(),
                     _CallersPassedTag(),
                     _CallersReportedListTag(),
                     _CallersPassedListTag(),
                     _SamplesReported(),
                     _SamplesPassed(),
                     _AlleleFreqAverageTag(),
                     _AlleleFreqRangeTag(),
                     _DepthAverageTag(),
                     _DepthRangeTag(),
                     _SomaticTag(),
                     _HCGenotypeTag()]

    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.add_tag_values(vcf_record)
        return vcf_record

    def get_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

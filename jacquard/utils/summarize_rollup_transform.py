"""Classes to summarize data for a sample-variant or variant as a whole.

A collection of individual Tag classes hold the metaheader and logic to
transform Jacquard-standardized VcfRecords.
"""
from __future__ import print_function, absolute_import, division

from decimal import Decimal
import re

import jacquard.utils.utils as utils
import jacquard.variant_caller_transforms.common_tags as common_tags


JQ_SUMMARY_TAG = "JQ_SUMMARY_"
JQ_REPORTED = "CALLERS_REPORTED_COUNT"
JQ_REPORTED_LIST = "CALLERS_REPORTED_LIST"
JQ_SAMPLES_REPORTED = "SAMPLES_REPORTED_COUNT"
JQ_PASSED = "CALLERS_PASSED_COUNT"
JQ_PASSED_LIST = "CALLERS_PASSED_LIST"
JQ_SAMPLES_PASSED = "SAMPLES_PASSED_COUNT"

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

def _get_non_null_values(record, sample, tag_name_regex):
    values = []
    try:
        for tag, value in list(record.sample_tag_values[sample].items()):
            if tag_name_regex.match(tag) and value != ".":
                values.append(value)
    except KeyError:
        raise utils.JQException("Sample [{}] was not recognized".format(sample))

    return values

def _average(numeric_values):
    average = str(sum(numeric_values)/Decimal(len(numeric_values)))
    return utils.round_two_digits(average)

def _range(numeric_values):
    if len(numeric_values) > 1:
        value_range = str(max(numeric_values) - min(numeric_values))
        return utils.round_two_digits(value_range)
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

def _add_caller_list_values(pattern, vcf_record, jq_global_variable):
    sample_tag = {}
    for sample, tags in vcf_record.sample_tag_values.items():
        sample_tag[sample] = []
        for tag in tags:
            pattern_match = pattern.match(tag)
            if pattern_match and tags[tag] != "0" and tags[tag] != ".":
                sample_tag[sample].append(pattern_match.group(1))

    for sample in sample_tag:
        sample_tag[sample] = ",".join(sorted(sample_tag[sample]))
        if not sample_tag[sample]:
            sample_tag[sample] = "."
    vcf_record.add_sample_tag_value(JQ_SUMMARY_TAG + jq_global_variable,
                                    sample_tag)

def _add_caller_count_values(pattern, vcf_record, jq_global_variable):
    sample_tag = {}
    for sample, tags in vcf_record.sample_tag_values.items():
        sample_tag[sample] = 0
        for tag in tags:
            if pattern.match(tag) and tags[tag] != "0" and tags[tag] != ".":
                sample_tag[sample] += int(tags[tag])
    vcf_record.add_sample_tag_value(JQ_SUMMARY_TAG + jq_global_variable,
                                    sample_tag)

def _add_sample_count_values(vcf_record,
                             jq_callers_global_variable,
                             jq_samples_global_variable):
    count = 0
    for tags in vcf_record.sample_tag_values.values():
        tag_key = JQ_SUMMARY_TAG + jq_callers_global_variable
        if tag_key in tags and tags[tag_key] != "0" and tags[tag_key] != ".":
            count += 1

    info_key = JQ_SUMMARY_TAG + jq_samples_global_variable
    vcf_record.add_info_field("=".join([info_key, str(count)]))

class _CallersReportedListTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.pattern = re.compile(r"JQ_(.*?)_{}"\
                                  .format(common_tags.CALLER_REPORTED_TAG))

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={}{},'
                'Number=.,'
                'Type=String,'
                #pylint: disable=line-too-long
                'Description="Comma-separated list variant callers which listed this variant in the Jacquard tagged VCF">')\
                .format(JQ_SUMMARY_TAG,
                        JQ_REPORTED_LIST)

    def add_tag_values(self, vcf_record):
        _add_caller_list_values(self.pattern, vcf_record, JQ_REPORTED_LIST)

class _CallersReportedTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.pattern = re.compile(r"JQ_(.*?)_{}"\
                                  .format(common_tags.CALLER_REPORTED_TAG))

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={}{},'
                'Number=1,'
                'Type=Integer,'
                #pylint: disable=line-too-long
                'Description="Count of variant callers which listed this variant in the Jacquard tagged VCF">')\
                .format(JQ_SUMMARY_TAG,
                        JQ_REPORTED)

    def add_tag_values(self, vcf_record):
        _add_caller_count_values(self.pattern, vcf_record, JQ_REPORTED)

class _CallersPassedListTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.pattern = re.compile(r"JQ_(.*?)_{}"\
                                  .format(common_tags.CALLER_PASSED_TAG))

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={}{},'
                'Number=.,'
                'Type=String,'
                #pylint: disable=line-too-long
                'Description="Comma-separated list of variant caller short-names where FILTER = PASS for this variant in the Jacquard tagged VCF">')\
                .format(JQ_SUMMARY_TAG,
                        JQ_PASSED_LIST)

    def add_tag_values(self, vcf_record):
        _add_caller_list_values(self.pattern, vcf_record, JQ_PASSED_LIST)

class _CallersPassedTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.pattern = re.compile(r"JQ_(.*?)_{}"\
                                  .format(common_tags.CALLER_PASSED_TAG))

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={}{},'
                'Number=1,'
                'Type=Integer,'
                #pylint: disable=line-too-long
                'Description="Count of variant callers where FILTER = PASS for this variant in the Jacquard tagged VCF">')\
                .format(JQ_SUMMARY_TAG,
                        JQ_PASSED)

    def add_tag_values(self, vcf_record):
        _add_caller_count_values(self.pattern, vcf_record, JQ_PASSED)

class _SamplesReported(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()

    @staticmethod
    def _get_metaheader():
        return ('##INFO=<ID={}{},'
                'Number=1,'
                'Type=Integer,'
                #pylint: disable=line-too-long
                'Description="Count of samples where this variant appeared in any of the Jacquard tagged VCFs (regardless of quality/filtering)">')\
                .format(JQ_SUMMARY_TAG,
                        JQ_SAMPLES_REPORTED)

    @staticmethod
    def add_tag_values(vcf_record):
        _add_sample_count_values(vcf_record, JQ_REPORTED, JQ_SAMPLES_REPORTED)

class _SamplesPassed(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()

    @staticmethod
    def _get_metaheader():
        return ('##INFO=<ID={}{},'
                'Number=1,'
                'Type=Integer,'
                #pylint: disable=line-too-long
                'Description="Count of samples where a variant caller passed the filter in any of the Jacquard tagged VCFs">')\
                .format(JQ_SUMMARY_TAG,
                        JQ_SAMPLES_PASSED)

    @staticmethod
    def add_tag_values(vcf_record):
        _add_sample_count_values(vcf_record, JQ_PASSED, JQ_SAMPLES_PASSED)

class _AlleleFreqRangeTag(object):
    #pylint: disable=too-few-public-methods
    _TAG_ID = "{}AF_RANGE".format(JQ_SUMMARY_TAG)
    _PATTERN = re.compile("^JQ_.*_AF$")

    def __init__(self):
        self.metaheader = self._get_metaheader()

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={},'
                'Number=1,'
                'Type=Float,'
                ##pylint: disable=line-too-long
                'Description="Max(allele frequency) - min (allele frequency) across recognized callers.">')\
                .format(_AlleleFreqRangeTag._TAG_ID)

    @staticmethod
    def add_tag_values(record):
        new_sample_tag_values = {}
        for sample in record.sample_tag_values:
            tag_values = _get_non_null_values(record,
                                              sample,
                                              _AlleleFreqRangeTag._PATTERN)

            aggregated_values = "."
            if tag_values:
                aggregated_values = _aggregate_numeric_values(tag_values,
                                                              _range)
            if not aggregated_values:
                #pylint: disable=line-too-long
                raise utils.JQException("Error summarizing values {} at record [{}:{} {}]",
                                        list(tag_values),
                                        record.chrom,
                                        record.pos,
                                        sample)

            new_sample_tag_values[sample] = aggregated_values

        record.add_sample_tag_value(_AlleleFreqRangeTag._TAG_ID,
                                    new_sample_tag_values)

class _AlleleFreqAverageTag(object):
    #pylint: disable=too-few-public-methods
    _TAG_ID = "{}AF_AVERAGE".format(JQ_SUMMARY_TAG)
    _PATTERN = re.compile("^JQ_.*_AF$")

    def __init__(self):
        self.metaheader = self._get_metaheader()

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={0},'
                'Number=1,'
                'Type=Float,'
                ##pylint: disable=line-too-long
                'Description="Average allele frequency across recognized variant callers that reported frequency for this position [average(JQ_*_AF)].">')\
                .format(_AlleleFreqAverageTag._TAG_ID)

    @staticmethod
    def add_tag_values(record):
        new_sample_tag_values = {}
        for sample in record.sample_tag_values:
            tag_values = _get_non_null_values(record,
                                              sample,
                                              _AlleleFreqAverageTag._PATTERN)

            aggregated_values = "."
            if tag_values:
                aggregated_values = _aggregate_numeric_values(tag_values,
                                                              _average)
            if not aggregated_values:
                #pylint: disable=line-too-long
                raise utils.JQException("Error summarizing values {} at record [{}:{} {}]",
                                        list(tag_values),
                                        record.chrom,
                                        record.pos,
                                        sample)

            new_sample_tag_values[sample] = aggregated_values

        record.add_sample_tag_value(_AlleleFreqAverageTag._TAG_ID,
                                    new_sample_tag_values)

class _DepthRangeTag(object):
    #pylint: disable=too-few-public-methods
    _TAG_ID = "{}DP_RANGE".format(JQ_SUMMARY_TAG)
    _PATTERN = re.compile("^JQ_.*_DP$")

    def __init__(self):
        self.metaheader = self._get_metaheader()

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={},'
                'Number=1,'
                'Type=Float,'
                ##pylint: disable=line-too-long
                'Description="Max(depth) - min (depth) across recognized callers.">'
                .format(_DepthRangeTag._TAG_ID))

    @staticmethod
    def add_tag_values(record):
        new_sample_tag_values = {}
        for sample in record.sample_tag_values:
            tag_values = _get_non_null_values(record,
                                              sample,
                                              _DepthRangeTag._PATTERN)

            aggregated_values = "."
            if tag_values:
                aggregated_values = _aggregate_numeric_values(tag_values,
                                                              _range)
            if not aggregated_values:
                #pylint: disable=line-too-long
                raise utils.JQException("Error summarizing values {} at record [{}:{} {}]",
                                        list(tag_values),
                                        record.chrom,
                                        record.pos,
                                        sample)

            new_sample_tag_values[sample] = aggregated_values

        record.add_sample_tag_value(_DepthRangeTag._TAG_ID,
                                    new_sample_tag_values)

class _DepthAverageTag(object):
    #pylint: disable=too-few-public-methods
    _TAG_ID = "{}DP_AVERAGE".format(JQ_SUMMARY_TAG)
    _PATTERN = re.compile("^JQ_.*_DP$")

    def __init__(self):
        self.metaheader = self._get_metaheader()

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={},'
                'Number=1,'
                'Type=Float,'
                ##pylint: disable=line-too-long
                'Description="Average allele frequency across recognized variant callers that reported frequency for this position; rounded to integer [round(average(JQ_*_DP))].">'
                .format(_DepthAverageTag._TAG_ID))

    @staticmethod
    def add_tag_values(record):
        new_sample_tag_values = {}
        for sample in record.sample_tag_values:
            tag_values = _get_non_null_values(record,
                                              sample,
                                              _DepthAverageTag._PATTERN)

            aggregated_values = "."
            if tag_values:
                aggregated_values = _aggregate_numeric_values(tag_values,
                                                              _average)
            if not aggregated_values:
                #pylint: disable=line-too-long
                raise utils.JQException("Error summarizing values {} at record [{}:{} {}]",
                                        list(tag_values),
                                        record.chrom,
                                        record.pos,
                                        sample)

            new_sample_tag_values[sample] = aggregated_values

        record.add_sample_tag_value(_DepthAverageTag._TAG_ID,
                                    new_sample_tag_values)

class _SomaticTag(object):
    #pylint: disable=too-few-public-methods
    _TAG_ID = "{}SOM_COUNT".format(JQ_SUMMARY_TAG)
    _PATTERN = re.compile("^JQ_.*_HC_SOM$")

    def __init__(self):
        self.metaheader = self._get_metaheader()

    @staticmethod
    def _get_metaheader():
        som_count = ('##FORMAT=<ID={},'
                     'Number=1,'
                     'Type=Integer,'
                     ##pylint: disable=line-too-long
                     'Description="Count of recognized variant callers that reported confident somatic call for this sample-position.">')\
                     .format(_SomaticTag._TAG_ID)
        return som_count

    @staticmethod
    def add_tag_values(record):
        new_sample_tag_values = {}
        for sample in record.sample_tag_values:
            tag_values = _get_non_null_values(record,
                                              sample,
                                              _SomaticTag._PATTERN)
            aggregated_values = "."
            if tag_values:
                aggregated_values = _aggregate_numeric_values(tag_values,
                                                              _count)
            if not aggregated_values:
                #pylint: disable=line-too-long
                raise utils.JQException("Error summarizing values {} at record [{}:{} {}]",
                                        list(tag_values),
                                        record.chrom,
                                        record.pos,
                                        sample)

            new_sample_tag_values[sample] = aggregated_values

        record.add_sample_tag_value(_SomaticTag._TAG_ID,
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
                     _SomaticTag()]

    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.add_tag_values(vcf_record)
        return vcf_record

    def get_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

from __future__ import print_function, absolute_import
import jacquard.utils as utils
import numpy as np
import re
from collections import defaultdict

JQ_CONSENSUS_TAG = "JQ_CONS_"
JQ_REPORTED = "CALLERS_REPORTED_COUNT"
JQ_REPORTED_LIST = "CALLERS_REPORTED_LIST"
JQ_PASSED = "CALLERS_PASSED_COUNT"
JQ_PASSED_LIST = "CALLERS_PASSED_LIST"

def _round_two_digits(value):
    split_value = value.split(".")

    if len(split_value[1]) <= 2:
        if split_value[1] == '0':
            return split_value[0]
        return value

    else:
        return str(round(100 * float(value))/100)

def _build_new_tags(vcf_record, tags, sample):
    desired_tags = []
    for tag in tags:
        desired_tag = vcf_record.sample_tag_values[sample][tag].split(",")

        #don't include null values in avg calculation
        altered_tag = [x for x in desired_tag if x != "."]
        if len(altered_tag) != 0:
            desired_tags.append(altered_tag)

    return desired_tags

def _get_tag_consensus_and_range(vcf_record, tags, all_ranges):
    tag_consensus = {}
    tag_range = {}

    for sample in vcf_record.sample_tag_values.keys():
        desired_tags = _build_new_tags(vcf_record, tags, sample)

        if len(desired_tags) == 0:
            tag_consensus[sample] = "."
        else:
            tag_consensus[sample] = _calculate_average(desired_tags)

        tag_range[sample] = _calculate_range(desired_tags, all_ranges)

    return tag_consensus, tag_range

def _calculate_average(tags):
    #pylint: disable=no-member
    tag_array = np.array(tags)
    rounded_tags = []

    for i in xrange(len(tag_array[0,])):
        tag_values = tag_array.astype(float)[:, i]
        rounded_tag = _round_two_digits(str(np.mean(tag_values)))
        rounded_tags.append(rounded_tag)

    return ",".join(rounded_tags)

def _calculate_range(tags, all_ranges):
    #don't calculate range if only called by one caller
    if len(tags) > 1:
        #pylint: disable=no-member
        cons_tag_array = np.array(tags)
        tag_range = []

        for i in xrange(len(cons_tag_array[0,])):
            tag_values = cons_tag_array.astype(float)[:, i]
            this_tag_range = np.max(tag_values) - np.min(tag_values)
            rounded_tag_range = _round_two_digits(str(this_tag_range))
            tag_range.append(rounded_tag_range)

            all_ranges.append(float(rounded_tag_range))
        return ",".join(tag_range)

    else:
        return "."

def _get_somatic_count(vcf_record, tags):
    somatic_count = {}

    for sample in vcf_record.sample_tag_values.keys():
        desired_tags = _build_new_tags(vcf_record, tags, sample)

        if len(desired_tags) == 0:
            somatic_count[sample] = "."
        else:
            somatic_count[sample] = _calculate_count(desired_tags)

    return somatic_count

def _calculate_count(tags):
    #pylint: disable=no-member
    tag_array = np.array(tags)
    counts = []

    for i in xrange(len(tag_array[0,])):
        tag_values = tag_array.astype(int)[:, i]
        count = sum(i == 1 for i in tag_values)
        counts.append(str(count))

    return ",".join(counts)

def _add_caller_list_values(pattern, vcf_record, jq_global_variable):
    sample_tag = defaultdict(list)
    for sample, tags in vcf_record.sample_tag_values.items():
        for tag in tags:
            pattern_match = pattern.match(tag)
            if pattern_match and tags[tag] != '0':
                sample_tag[sample].append(pattern_match.group(1))

    for sample in sample_tag:
        sample_tag[sample] = ",".join(sorted(sample_tag[sample]))
    vcf_record.add_sample_tag_value(JQ_CONSENSUS_TAG + jq_global_variable,
                                    sample_tag)

def _add_caller_count_values(pattern, vcf_record, jq_global_variable):
    sample_tag = {}
    for sample, tags in vcf_record.sample_tag_values.items():
        sample_tag[sample] = 0
        for tag in tags:
            if pattern.match(tag):
                sample_tag[sample] += int(tags[tag])
    vcf_record.add_sample_tag_value(JQ_CONSENSUS_TAG + jq_global_variable,
                                    sample_tag)

class _CallersReportedListTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.all_ranges = []
        self.name = ""
        self.pattern = re.compile(r"JQ_(.*)_{}".format(JQ_REPORTED))

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={}{},'\
                'Number=.,'\
                'Type=String,'\
                'Description="Comma-separated list variant callers '\
                'which listed this variant in the Jacquard tagged VCF",'\
                'Source="Jacquard",'\
                'Version="{}">'\
                .format(JQ_CONSENSUS_TAG,
                        JQ_REPORTED_LIST,
                        utils.__version__))

    def add_tag_values(self, vcf_record):
        _add_caller_list_values(self.pattern, vcf_record, JQ_REPORTED_LIST)

class _CallersReportedTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.all_ranges = []
        self.name = ""
        self.pattern = re.compile(r"JQ_(.*)_{}".format(JQ_REPORTED))

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={}{},'\
                'Number=1,'\
                'Type=Integer,'\
                'Description="Count of variant callers which listed '\
                'this variant in the Jacquard tagged VCF",'\
                'Source="Jacquard",'\
                'Version="{}">\n')\
                .format(JQ_CONSENSUS_TAG,
                        JQ_REPORTED,
                        utils.__version__)

    def add_tag_values(self, vcf_record):
        _add_caller_count_values(self.pattern, vcf_record, JQ_REPORTED)

class _CallersPassedListTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.all_ranges = []
        self.name = ""
        self.pattern = re.compile(r"JQ_(.*)_{}".format(JQ_PASSED))

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={}{},'\
                'Number=.,'\
                'Type=String,'\
                'Description="Comma-separated list of variant caller '\
                'short-names where FILTER = PASS for this variant in '\
                'the Jacquard tagged VCF",'\
                'Source="Jacquard",'\
                'Version="{}">'\
                .format(JQ_CONSENSUS_TAG,
                        JQ_PASSED_LIST,
                        utils.__version__))

    def add_tag_values(self, vcf_record):
        _add_caller_list_values(self.pattern, vcf_record, JQ_PASSED_LIST)

class _CallersPassedTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.all_ranges = []
        self.name = ""
        self.pattern = re.compile(r"JQ_(.*)_{}".format(JQ_PASSED))

    @staticmethod
    def _get_metaheader():
        return ('##FORMAT=<ID={}{},'\
                'Number=1,Type=Integer,'\
                'Description="Count of variant callers where FILTER = PASS '\
                'for this variant in the Jacquard tagged VCF",'\
                'Source="Jacquard",'\
                'Version="{}">\n')\
                .format(JQ_CONSENSUS_TAG,
                        JQ_PASSED,
                        utils.__version__)

    def add_tag_values(self, vcf_record):
        _add_caller_count_values(self.pattern, vcf_record, JQ_PASSED)

class _AlleleFreqTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.all_ranges = []
        self.name = "AF"

    @staticmethod
    def _get_metaheader():
        af_average = ('##FORMAT=<ID={0}AF_AVERAGE,'
                      'Number=1,'
                      'Type=Float,'
                      ##pylint: disable=line-too-long
                      'Description="Average allele frequency across recognized variant callers that reported frequency for this position [average(JQ_*_AF)].",'
                      'Source="Jacquard",'
                      'Version="{1}">').format(JQ_CONSENSUS_TAG,
                                               utils.__version__)
        af_range = ('##FORMAT=<ID={0}AF_RANGE,'
                    'Number=1,'
                    'Type=Float,'
                    ##pylint: disable=line-too-long
                    'Description="Max(allele frequency) - min (allele frequency) across recognized callers.",'
                    'Source="Jacquard",'
                    'Version="{1}>">').format(JQ_CONSENSUS_TAG,
                                              utils.__version__)
        return "\n".join([af_average, af_range])

    @staticmethod
    def _get_allele_freq_tags(vcf_record):
        tags = []
        for tag in vcf_record.format_tags:
            if tag.startswith("JQ_") and tag.endswith("_AF"):
                tags.append(tag)

        return tags

    def add_tag_values(self, vcf_record):
        tags = self._get_allele_freq_tags(vcf_record)
        tag_consensus, tag_range = _get_tag_consensus_and_range(vcf_record,
                                                                tags,
                                                                self.all_ranges)

        vcf_record.add_sample_tag_value(JQ_CONSENSUS_TAG + "AF_AVERAGE",
                                        tag_consensus)
        vcf_record.add_sample_tag_value(JQ_CONSENSUS_TAG + "AF_RANGE",
                                        tag_range)


class _DepthTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.all_ranges = []
        self.name = "DP"

    @staticmethod
    def _get_metaheader():
        dp_average = ('##FORMAT=<ID={0}DP_AVERAGE,'
                      'Number=1,'
                      'Type=Float,'
                      ##pylint: disable=line-too-long
                      'Description="Average allele frequency across recognized variant callers that reported frequency for this position; rounded to integer [round(average(JQ_*_DP))].",'
                      'Source="Jacquard",'
                      'Version="{1}">').format(JQ_CONSENSUS_TAG,
                                               utils.__version__)
        dp_range = ('##FORMAT=<ID={0}DP_RANGE,'
                    'Number=1,'
                    'Type=Float,'
                    ##pylint: disable=line-too-long
                    'Description="Max(depth) - min (depth) across recognized callers.",'
                    'Source="Jacquard",'
                    'Version="{1}>">').format(JQ_CONSENSUS_TAG,
                                              utils.__version__)
        return "\n".join([dp_average, dp_range])

    @staticmethod
    def _get_depth_tags(vcf_record):
        tags = []
        for tag in vcf_record.format_tags:
            if tag.startswith("JQ_") and tag.endswith("_DP"):
                tags.append(tag)

        return tags

    def add_tag_values(self, vcf_record):
        tags = self._get_depth_tags(vcf_record)
        tag_consensus, tag_range = _get_tag_consensus_and_range(vcf_record,
                                                                tags,
                                                                self.all_ranges)

        vcf_record.add_sample_tag_value(JQ_CONSENSUS_TAG + "DP_AVERAGE",
                                        tag_consensus)
        vcf_record.add_sample_tag_value(JQ_CONSENSUS_TAG + "DP_RANGE",
                                        tag_range)


class _SomaticTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.all_ranges = []
        self.name = "SOM"

    @staticmethod
    def _get_metaheader():
        som_count = ('##FORMAT=<ID={0}SOM_COUNT,'
                     'Number=1,'
                     'Type=Integer,'
                     ##pylint: disable=line-too-long
                     'Description="Count of recognized variant callers that reported confident somatic call for this sample-position.",'
                     'Source="Jacquard",'
                     'Version="{1}">').format(JQ_CONSENSUS_TAG,
                                              utils.__version__)
        return som_count

    @staticmethod
    def _get_somatic_tags(vcf_record):
        tags = []
        for tag in vcf_record.format_tags:
            if tag.startswith("JQ_") and tag.endswith("_HC_SOM"):
                tags.append(tag)

        return tags

    def add_tag_values(self, vcf_record):
        tags = self._get_somatic_tags(vcf_record)
        somatic_count = _get_somatic_count(vcf_record, tags)

        vcf_record.add_sample_tag_value(JQ_CONSENSUS_TAG + "SOM_COUNT",
                                        somatic_count)


class ConsensusCaller(object):
    def __init__(self):
        self.tags = [_CallersReportedTag(),
                     _CallersPassedTag(),
                     _AlleleFreqTag(),
                     _DepthTag(),
                     _SomaticTag()]
        self.ranges = {}

    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.add_tag_values(vcf_record)
            self.ranges[tag.name] = tag.all_ranges

        return vcf_record.asText()

    def get_consensus_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

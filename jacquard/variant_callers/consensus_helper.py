# pylint: disable=C0111
import numpy as np

import jacquard.utils as utils

JQ_CONSENSUS_TAG = "JQ_CONS_"

def _round_two_digits(value):
    split_value = value.split(".")

    if len(split_value[1]) <= 2:
        if split_value[1] == '0':
            return split_value[0]
        return value

    else:
        return str(round(100 * float(value))/100)

def _build_desired_tags(vcf_record, tags, sample):
    desired_tags = []
    for tag in tags:
        desired_tag = vcf_record.sample_dict[sample][tag].split(",")

        #don't include null values in avg calculation
        altered_tag = [x for x in desired_tag if x != "."]
        if len(altered_tag) != 0:
            desired_tags.append(altered_tag)

    return desired_tags

def _get_tag_consensus_and_range(vcf_record, tags, all_ranges):
    tag_consensus = {}
    tag_range = {}

    for sample in vcf_record.sample_dict.keys():
        desired_tags = _build_desired_tags(vcf_record, tags, sample)

        if len(desired_tags) == 0:
            tag_consensus[sample] = "."
        else:
            tag_consensus[sample] = _calculate_average(desired_tags)

        tag_range[sample] = _calculate_range(desired_tags, all_ranges)

    return tag_consensus, tag_range

def _calculate_average(tags):
    tag_array = np.array(tags)
    rounded_tags = []

    for i in xrange(len(tag_array[0,])):
        tag_values = tag_array.astype(float)[:,i]
        rounded_tag = _round_two_digits(str(np.mean(tag_values)))
        rounded_tags.append(rounded_tag)

    return ",".join(rounded_tags)

def _calculate_range(tags, all_ranges):
    #don't calculate range if only called by one caller
    if len(tags) > 1:
        cons_tag_array = np.array(tags)
        tag_range = []

        for i in xrange(len(cons_tag_array[0,])):
            tag_values = cons_tag_array.astype(float)[:,i]
            this_tag_range = np.max(tag_values) - np.min(tag_values)
            rounded_tag_range = _round_two_digits(str(this_tag_range))
            tag_range.append(rounded_tag_range)

            all_ranges.append(float(rounded_tag_range))
        return ",".join(tag_range)

    else:
        return "."

def _calculate_population_values(all_ranges):
    for ranges in all_ranges.values():
        pop_mean_range = str(sum(ranges)/len(ranges))
        pop_std_range = str(np.std(ranges))

        rounded_pop_mean_range = _round_two_digits(pop_mean_range)
        rounded_pop_std_range = _round_two_digits(pop_std_range)

        pop_mean_range = float(rounded_pop_mean_range)
        pop_std_range = float(rounded_pop_std_range)

        return (pop_mean_range, pop_std_range)

def _calculate_zscore(vcf_record, tag, pop_mean_range, pop_std_range):
    zscore_dict = {}
    for sample in vcf_record.sample_dict.keys():
        samp_range = vcf_record.sample_dict[sample][tag]
        if samp_range != ".":
            samp_range = float(samp_range)
            if pop_std_range != 0.0:
                zscore = (samp_range - pop_mean_range)/pop_std_range
            else:
                zscore = "."
            zscore_dict[sample] = zscore
        else:
            zscore_dict[sample] = "."

    return zscore_dict

def _get_somatic_count(vcf_record, tags):
    somatic_count = {}

    for sample in vcf_record.sample_dict.keys():
        desired_tags = _build_desired_tags(vcf_record, tags, sample)

        if len(desired_tags) == 0:
            somatic_count[sample] = "."
        else:
            somatic_count[sample] = _calculate_count(desired_tags)

    return somatic_count

def _calculate_count(tags):
    tag_array = np.array(tags)
    counts = []

    for i in xrange(len(tag_array[0,])):
        tag_values = tag_array.astype(int)[:,i]
        count = sum(i == 1 for i in tag_values)
        counts.append(str(count))

    return ",".join(counts)

class _AlleleFreqTag():
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.all_ranges = []
        self.name = "AF"

    def _get_metaheader(self):
        af_average = '##FORMAT=<ID={0}AF_AVERAGE,Number=1,Type=Float,' \
                      'Description="Average allele frequency across ' \
                      'recognized variant callers that reported ' \
                      'frequency for this position [average(JQ_*_AF)].",' \
                      'Source="Jacquard",Version="{1}">'\
                      .format(JQ_CONSENSUS_TAG, utils.__version__)
        af_range = '##FORMAT=<ID={0}AF_RANGE, Number=1,Type=Float,' \
                   'Description="Max(allele frequency) - min (allele '\
                   'frequency) across recognized callers.",Source="Jacquard",'\
                   'Version="<{1}>">'\
                   .format(JQ_CONSENSUS_TAG, utils.__version__)
        af_zscore = '##FORMAT=<ID={0}AF_ZSCORE,Number=1,Type=Float,'\
                    'Description="Jacquard measure of concordance of reported '\
                    'allele frequencies across callers. [(this AF range - '\
                    'mean AF range)/standard dev(all AF ranges)]. If '\
                    'consensus value from <2 values will be [.]",Source="'\
                    'Jacquard",Version="<{1}>"'\
                    .format(JQ_CONSENSUS_TAG, utils.__version__)
        return "\n".join([af_average, af_range, af_zscore])

    def _get_allele_freq_tags(self, vcf_record):
        tags = []
        for tag in vcf_record.format_set:
            if tag.startswith("JQ_") and tag.endswith("_AF"):
                tags.append(tag)

        return tags

    def insert_consensus(self, vcf_record):
        tags = self._get_allele_freq_tags(vcf_record)
        tag_consensus, tag_range = _get_tag_consensus_and_range(vcf_record,
                                                                tags,
                                                                self.all_ranges)

        vcf_record.insert_format_field(JQ_CONSENSUS_TAG + "AF_AVERAGE",
                                      tag_consensus)
        vcf_record.insert_format_field(JQ_CONSENSUS_TAG + "AF_RANGE",
                                      tag_range)

    def insert_zscore(self, vcf_record, pop_mean_range, pop_std_range):
        tag = JQ_CONSENSUS_TAG + "AF_RANGE"
        zscore_dict = _calculate_zscore(vcf_record,
                                       tag,
                                       pop_mean_range,
                                       pop_std_range)

        vcf_record.insert_format_field(JQ_CONSENSUS_TAG + "AF_ZSCORE",
                                       zscore_dict)

class _DepthTag():
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.all_ranges = []
        self.name = "DP"

    def _get_metaheader(self):
        dp_average = '##FORMAT=<ID={0}DP_AVERAGE,Number=1,Type=Float,' \
                      'Description="Average allele frequency across ' \
                      'recognized variant callers that reported ' \
                      'frequency for this position; rounded to integer ' \
                      '[round(average(JQ_*_DP))].",' \
                      'Source="Jacquard",Version="{1}">'\
                      .format(JQ_CONSENSUS_TAG, utils.__version__)
        dp_range = '##FORMAT=<ID={0}DP_RANGE, Number=1,Type=Float,' \
                   'Description="Max(depth) - min (depth) '\
                   'across recognized callers.",Source="Jacquard",'\
                   'Version="<{1}>">'\
                   .format(JQ_CONSENSUS_TAG, utils.__version__)
        dp_zscore = '##FORMAT=<ID={0}DP_ZSCORE,Number=1,Type=Float,'\
                    'Description="Jacquard measure of concordance of reported '\
                    'depths across callers. [(this DP range - '\
                    'mean DP range)/standard dev(all DP ranges)]. If '\
                    'consensus value from <2 values will be [.]",Source="'\
                    'Jacquard",Version="<{1}>"'\
                    .format(JQ_CONSENSUS_TAG, utils.__version__)
        return "\n".join([dp_average, dp_range, dp_zscore])

    def _get_depth_tags(self, vcf_record):
        tags = []
        for tag in vcf_record.format_set:
            if tag.startswith("JQ_") and tag.endswith("_DP"):
                tags.append(tag)

        return tags

    def insert_consensus(self, vcf_record):
        tags = self._get_depth_tags(vcf_record)
        tag_consensus, tag_range = _get_tag_consensus_and_range(vcf_record,
                                                                tags,
                                                                self.all_ranges)

        vcf_record.insert_format_field(JQ_CONSENSUS_TAG + "DP_AVERAGE",
                                      tag_consensus)
        vcf_record.insert_format_field(JQ_CONSENSUS_TAG + "DP_RANGE",
                                       tag_range)

    def insert_zscore(self, vcf_record, pop_mean_range, pop_std_range):
        tag = JQ_CONSENSUS_TAG + "DP_RANGE"
        zscore_dict = _calculate_zscore(vcf_record,
                                       tag,
                                       pop_mean_range,
                                       pop_std_range)

        vcf_record.insert_format_field(JQ_CONSENSUS_TAG + "DP_ZSCORE",
                                       zscore_dict)

class _SomaticTag():
    def __init__(self):
        self.metaheader = self._get_metaheader()
        self.all_ranges = []
        self.name = "SOM"

    def _get_metaheader(self):
        som_count = '##FORMAT=<ID={0}SOM_COUNT,Number=1,Type=Integer,' \
                      'Description="Count of recognized variant callers,' \
                      'which reported confident somatic call for this'\
                      'sample-position.",Source="Jacquard",Version="{1}">'\
                      .format(JQ_CONSENSUS_TAG, utils.__version__)
        return som_count

    def _get_somatic_tags(self, vcf_record):
        tags = []
        for tag in vcf_record.format_set:
            if tag.startswith("JQ_") and tag.endswith("_HC_SOM"):
                tags.append(tag)

        return tags

    def insert_consensus(self, vcf_record):
        tags = self._get_somatic_tags(vcf_record)
        somatic_count = _get_somatic_count(vcf_record, tags)

        vcf_record.insert_format_field(JQ_CONSENSUS_TAG + "SOM_COUNT",
                                      somatic_count)

class ConsensusHelper():
    def __init__(self):
        self.tags = [_AlleleFreqTag(), _DepthTag()]
        self.ranges = {}

    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.insert_consensus(vcf_record)
            self.ranges[tag.name] = tag.all_ranges

        return vcf_record.asText()

    def get_population_values(self):
        pop_values = {}
        for tag in self.tags:
            (pop_mean_range,
             pop_std_range) = _calculate_population_values(self.ranges)
            pop_values[tag.name] = [pop_mean_range, pop_std_range]

        return pop_values

    def add_zscore(self, vcf_record, pop_values):
        for tag in self.tags:
            pop_mean_range, pop_std_range = pop_values[tag.name]
            tag.insert_zscore(vcf_record, pop_mean_range, pop_std_range)

        return vcf_record.asText()

    def get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

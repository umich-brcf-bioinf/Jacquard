# pylint: disable=C0111
from collections import defaultdict
import numpy as np

import jacquard.utils as utils

JQ_CONSENSUS_TAG = "JQ_CONS_"

class _AlleleFreqTag():
    def __init__(self):
        self.metaheader = self.get_metaheader()
        self.all_ranges = []
        self.name = "af"

    def get_metaheader(self):
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

    def format(self, vcf_record):
        cons_freqs = {}
        range_freqs = {}
        tags = self._get_allele_freq_tags(vcf_record)

        for sample in vcf_record.sample_dict.keys():
            freqs = []
            for tag in tags:
                freq = vcf_record.sample_dict[sample][tag].split(",")

                #don't include null values in avg calculation
                altered_freq = [x for x in freq if x != "."]
                if len(altered_freq) != 0:
                    freqs.append(altered_freq)

            if len(freqs) == 0:
                cons_freqs[sample] = "."
            else:
                self._validate_multAlts(freqs)
                avg_freqs = self._calculate_average(freqs)
                cons_freqs[sample] = avg_freqs

            range_freqs[sample] = self._calculate_range(freqs)

        vcf_record.insert_format_field(JQ_CONSENSUS_TAG + "AF_AVERAGE",
                                      cons_freqs)
        vcf_record.insert_format_field(JQ_CONSENSUS_TAG + "AF_RANGE",
                                      range_freqs)

    def _validate_multAlts(self, freqs):
        length = len(freqs[0])
        for freq in freqs:
            if len(freq) != length:
                raise utils.JQException("Inconsistent number of mult-alts "
                                        "found in VCF file.")

    def _get_allele_freq_tags(self, vcfRecord):
        tags = []
        for tag in vcfRecord.format_set:
            if tag.startswith("JQ_") and tag.endswith("_AF"):
                tags.append(tag)

        return tags

    def _calculate_average(self, freqs):
        freq_array = np.array(freqs)
        rounded_freqs = []

        for i in xrange(len(freq_array[0,])):
            freq_values = freq_array.astype(float)[:,i]
            rounded_freq = self._roundTwoDigits(str(np.mean(freq_values)))
            rounded_freqs.append(rounded_freq)

        return ",".join(rounded_freqs)

    def _roundTwoDigits(self, value):
        split_value = value.split(".")

        if len(split_value[1]) <= 2:
            if split_value[1] == '0':
                return split_value[0]
            return value

        else:
            return str(round(100 * float(value))/100)

    def _calculate_range(self, freqs):
        #don't calculate range if only called by one caller
        if len(freqs) > 1:
            cons_freq_array = np.array(freqs)
            af_range = []

            for i in xrange(len(cons_freq_array[0,])):
                freq_values = cons_freq_array.astype(float)[:,i]
                this_af_range = np.max(freq_values) - np.min(freq_values)
                rounded_af_range = self._roundTwoDigits(str(this_af_range))
                af_range.append(rounded_af_range)

                self.all_ranges.append(float(rounded_af_range))
            return ",".join(af_range)

        else:
            return "."

    def calculate_pop_values(self, all_ranges):
        for ranges in all_ranges.values():
            pop_mean_range = str(sum(ranges)/len(ranges))
            pop_std_range = str(np.std(ranges))

            rounded_pop_mean_range = self._roundTwoDigits(pop_mean_range)
            rounded_pop_std_range = self._roundTwoDigits(pop_std_range)

            pop_mean_range = float(rounded_pop_mean_range)
            pop_std_range = float(rounded_pop_std_range)

            return (pop_mean_range, pop_std_range)

    def calculate_zscore(self, vcf_record, pop_mean_range, pop_std_range):
        tag = JQ_CONSENSUS_TAG + "AF_RANGE"
        zscore_dict = {}

        for sample in vcf_record.sample_dict.keys():
            samp_range = vcf_record.sample_dict[sample][tag]
            if samp_range != ".":
                samp_range = float(samp_range)
                zscore = (samp_range - pop_mean_range)/pop_std_range if pop_std_range != 0.0 else "."
                zscore_dict[sample] = zscore
            else:
                zscore_dict[sample] = "."

        vcf_record.insert_format_field(JQ_CONSENSUS_TAG + "AF_ZSCORE",
                                       zscore_dict)

class ConsensusHelper():
    def __init__(self):
        self.tags = [_AlleleFreqTag()]
        self.ranges = {}

    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.format(vcf_record)
            self.ranges[tag.name] = tag.all_ranges

        return vcf_record.asText()

    def get_population_values(self):
        pop_values ={}
        for tag in self.tags:
            (pop_mean_range, pop_std_range) = tag.calculate_pop_values(self.ranges)
            pop_values[tag.name] = [pop_mean_range, pop_std_range]

        return pop_values

    def add_zscore(self, vcf_record, pop_values):
        for tag in self.tags:
            pop_mean_range, pop_std_range = pop_values[tag.name]
            tag.calculate_zscore(vcf_record, pop_mean_range, pop_std_range)

        return vcf_record.asText()

    def get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]


# pylint: disable=C0111
import numpy as np

import jacquard.utils as utils

JQ_CONSENSUS_TAG = "JQ_CONS_"

class _AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}AF_AVERAGE,Number=1,Type=Float,' \
                          'Description="Average allele frequency across ' \
                          'recognized variant callers that reported ' \
                          'frequency for this position [average(JQ_*_AF)].", ' \
                          'Source="Jacquard", Version="{1}">'\
                          .format(JQ_CONSENSUS_TAG, utils.__version__)

    def format(self, vcfRecord):
        cons_freqs = {}
        tags = self._get_allele_freq_tags(vcfRecord)

        for sample in vcfRecord.sample_dict.keys():
            freqs = []
            for tag in tags:
                freq = vcfRecord.sample_dict[sample][tag].split(",")
                freqs.append(freq)

            if len(freqs) == 0:
                cons_freqs[sample] = "."
            else:
                self._validate_multAlts(freqs)
    
                freqs = self._calculate_average(freqs)
                cons_freqs[sample] = freqs

        vcfRecord.insert_format_field(JQ_CONSENSUS_TAG + "AF_AVERAGE",
                                      cons_freqs)

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


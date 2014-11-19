'''
Created on Nov 18, 2014

@author: kmeng
'''

import jacquard.utils as utils
import re
import os
import numpy as np
from collections import defaultdict

JQ_CONSENSUS_TAG = "JQ_CONS_"

class _AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}AF_AVERAGE,Number=1,Type=Float,Description="Average allele frequency across recognized variant callers that reported frequency for this position [average(JQ_*_AF)].", Source="Jacquard", Version="{1}">'.format(JQ_CONSENSUS_TAG, utils.__version__)

    def format(self, vcfRecord):
        cons_freqs = {}
        tags = []
        for tag in vcfRecord.format_set:
            if tag.startswith("JQ_") and tag.endswith("_AF"):
                tags.append(tag)
        for key in vcfRecord.sample_dict.keys():
            freqs = []
            for tag in tags:
                freq = vcfRecord.sample_dict[key][tag].split(",")
                freqs.append(freq)
            freqs = self._averager(freqs)
            cons_freqs[key] = freqs
        vcfRecord.insert_format_field(JQ_CONSENSUS_TAG + "AF_AVERAGE",cons_freqs)

    def _averager(self, freqs):
        freq_array = np.array(freqs)
        if freq_array.dtype != "S1": # Hacky way to check to see if the numpy arrays are all aligned.
            raise utils.JQException("Inconsistent number of multiple Alts found in VCF file.")
        new_freqs = []
        for i in xrange(len(freq_array[0,])):
            new_freq = self._roundTwoDigits(str(np.mean(freq_array.astype(float)[:,i])))
            new_freqs.append(new_freq)
        return ",".join(new_freqs)
    
    def _roundTwoDigits(self, value): 
        split_value = value.split(".")
        if len(split_value[1]) <= 2:
            if split_value[1] == '0':
                return split_value[0]
            return value
        else:
            return str(round(100 * float(value))/100)


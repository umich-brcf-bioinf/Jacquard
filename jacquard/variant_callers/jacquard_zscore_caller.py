#pylint: disable=too-few-public-methods, unused-argument
from __future__ import print_function, absolute_import
import math
from jacquard import __version__

_JQ_CONSENSUS_TAG = "JQ_CONS_"

#TODO: (cgates): Define tag ids as public, class-level constants so dependent tags can reference them directly
class _AlleleFreqZScoreTag(object):
    TAG_ID = "{0}AF_ZSCORE".format(_JQ_CONSENSUS_TAG)
    _RANGE_TAG = "{0}AF_RANGE".format(_JQ_CONSENSUS_TAG)
    _METAHEADER_DESCRIPTION = ('''
Concordance of reported allele frequencies across callers: 
[(this AF range - mean AF range)/standard dev(all AF ranges)]. 
Values with null or missing AF range will be assigned zscore of \'.\'; 
for multi-valued ranges, zscore is of largest range.''').replace("\n", "")

    def __init__(self, vcf_reader):
        self.tag_id = self.TAG_ID
        self.tag = _ZScoreTag(self.TAG_ID,
                              self._METAHEADER_DESCRIPTION,
                              self._RANGE_TAG,
                              vcf_reader)

    @property
    def metaheaders(self):
        return self.tag.metaheaders

    def add_tag_values(self, vcf_record):
        self.tag.add_tag_values(vcf_record)

class _DepthZScoreTag(object):
    TAG_ID = "{0}DP_ZSCORE".format(_JQ_CONSENSUS_TAG)
    _RANGE_TAG = "{0}DP_RANGE".format(_JQ_CONSENSUS_TAG)
    _METAHEADER_DESCRIPTION = ('''
Concordance of reported depth across callers: 
[(this DP range - mean DP range)/standard dev(all DP ranges)]. 
Values with null or missing DP range will be assigned zscore 
of \'.\'.''').replace("\n", "")

    def __init__(self, vcf_reader):
        self.tag_id = self.TAG_ID
        self.tag = _ZScoreTag(self.TAG_ID,
                              self._METAHEADER_DESCRIPTION,
                              self._RANGE_TAG,
                              vcf_reader)

    @property
    def metaheaders(self):
        return self.tag.metaheaders

    def add_tag_values(self, vcf_record):
        self.tag.add_tag_values(vcf_record)

class _ZScoreTag(object):
    _STRING_FORMAT = "{0:.2f}"
    _EXECUTION_FORMAT = "##jacquard.consensus.{0}.{1}_{2}={3}"
    _METAHEADER_FORMAT = ('##FORMAT=<ID={0},'
                          'Number=1,'
                          'Type=Float,'
                          'Description="{1}",'
                          'Source="Jacquard",'
                          'Version="{2}">')

    def __init__(self,
                 tag_id,
                 metaheader_description,
                 dependent_tag_id,
                 vcf_reader):
        self._tag_id = tag_id
        self._dependent_tag_id = dependent_tag_id
        self._mean, self._stdev = self._init_population_stats(vcf_reader,
                                                              dependent_tag_id)
        self._metaheaders = self._init_metaheaders(tag_id,
                                                   metaheader_description,
                                                   dependent_tag_id,
                                                   self._mean,
                                                   self._stdev)

    def _init_metaheaders(self,
                          tag_id,
                          metaheader_description,
                          dependent_tag_id,
                          mean,
                          stdev):
        #pylint: disable=too-many-arguments
        metaheaders = []
        metaheaders.append(self._EXECUTION_FORMAT.format(tag_id,
                                                         dependent_tag_id,
                                                         "mean",
                                                         mean))
        metaheaders.append(self._EXECUTION_FORMAT.format(tag_id,
                                                         dependent_tag_id,
                                                         "stdev",
                                                         stdev))
        tag_metaheader = self._METAHEADER_FORMAT.format(tag_id,
                                                        metaheader_description,
                                                        __version__)
        metaheaders.append(tag_metaheader)
        return tuple(metaheaders)

    @property
    def metaheaders(self):
        return self._metaheaders


    def _ok_to_add_tag_values(self, vcf_record):
        return self._stdev and self._dependent_tag_id in vcf_record.format_tags


    def _zscore_as_str(self, zscore):
        if zscore == ".":
            return zscore
        else:
            return self._STRING_FORMAT.format(zscore)

    def add_tag_values(self, vcf_record):
        if not self._ok_to_add_tag_values(vcf_record):
            return
        sample_values = {}
        for sample_name in vcf_record.sample_tag_values:
            zscore = "."
            tag_values = vcf_record.sample_tag_values[sample_name]
            value = self._get_dependent_value(tag_values,
                                              self._dependent_tag_id)
            if  value is not None:
                zscore = (value - self._mean) / self._stdev
            sample_values[sample_name] = self._zscore_as_str(zscore)
        vcf_record.add_sample_tag_value(self._tag_id,
                                        sample_values)

    @staticmethod
    def _get_dependent_value(tag_values, dependent_tag_id):
        try:
            values = tag_values[dependent_tag_id].split(",")
            return max([float(value) for value in values])
        except KeyError:
            return None
        except ValueError:
            return None

    def _init_population_stats(self, vcf_reader, dependent_tag_id):
        #pylint: disable=invalid-name
        n = 0
        mean = 0
        M2 = 0
        try:
            vcf_reader.open()
            for vcf_record in vcf_reader.vcf_records():
                for tag_values in vcf_record.sample_tag_values.values():
                    value = self._get_dependent_value(tag_values,
                                                      dependent_tag_id)
                    if value is not None:
                        n += 1
                        delta = value - mean
                        mean += delta/n
                        M2 += delta * (value - mean)
        finally:
            vcf_reader.close()

        stdev = 0
        if n == 0:
            mean = None
            stdev = None
        elif n >= 2:
            variance = M2/n
            stdev = math.sqrt(variance)

        return mean, stdev


class ZScoreCaller(object):
    def __init__(self, vcf_reader):
        self._tags = [_AlleleFreqZScoreTag(vcf_reader),
                      _DepthZScoreTag(vcf_reader)]
        self._metaheaders = self._init_metaheaders(self._tags)

    @staticmethod
    def _init_metaheaders(tags):
        metaheaders = []
        for tag in tags:
            metaheaders.extend(tag.metaheaders)
        return tuple(metaheaders)

    @property
    def metaheaders(self):
        return self._metaheaders

    def add_tags(self, vcf_record):
        for tag in self._tags:
            tag.add_tag_values(vcf_record)
        return vcf_record.asText()


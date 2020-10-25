"""Interprets MuTect VCF files adding Jacquard standard information.

MuTect VCFs are assumed to have a ".vcf" extension and have a valid
"##MuTect=..." metaheader.
"""
from __future__ import print_function, absolute_import, division

import argparse
import re

import jacquard.utils.utils as utils
import jacquard.variant_caller_transforms.common_tags as common_tags
import jacquard.utils.vcf as vcf

JQ_MUTECT_TAG = "JQ_MT_"
MUTECT_ABBREVIATION = "MT"
VERSION = "v1.1 - v2.2"

class _GenotypeTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods
    @classmethod
    def _standardize_gt(cls, value):
        if value == "0":
            value = "0/0"
        return value

    def __init__(self):
        super(self.__class__,
              self).__init__(MUTECT_ABBREVIATION,
                             common_tags.GENOTYPE_TAG,
                             "Jacquard genotype (based on GT)")

    def add_tag_values(self, vcf_record):
        if "GT" in vcf_record.format_tags:
            sample_values = {}
            for samp in vcf_record.sample_tag_values:
                genotype = vcf_record.sample_tag_values[samp]["GT"]
                sample_values[samp] = _GenotypeTag._standardize_gt(genotype)
            vcf_record.add_sample_tag_value(self.tag_id, sample_values)

class _AlleleFreqTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods
    def __init__(self, metaheaders):
        self.source_tag = self._determine_source_tag(metaheaders)
        super(self.__class__,
              self).__init__(MUTECT_ABBREVIATION,
                             common_tags.ALLELE_FREQ_TAG,
                             ('Jacquard allele frequency for MuTect: '
                              'Decimal allele frequency rounded to 4 digits '
                              '(based on {})').format(self.source_tag))

    @staticmethod
    def _determine_source_tag(metaheaders):
        for header in metaheaders:
            if header.startswith('##FORMAT=<ID=AF,'):
                return 'AF'
            if header.startswith('##FORMAT=<ID=FA,'):
                return 'FA'
        msg = ('could not determine the correct allele frequency '
               'FORMAT tag in the source MuTect file')
        raise utils.JQException(msg)

    def add_tag_values(self, vcf_record):
        if self.source_tag in vcf_record.format_tags:
            sample_values = {}
            for sample in vcf_record.sample_tag_values:
                freq = vcf_record.sample_tag_values[sample][self.source_tag].split(",")
                sample_values[sample] = self._standardize_af(freq)
            vcf_record.add_sample_tag_value(self.tag_id, sample_values)

    @staticmethod
    def _standardize_af(value):
        new_values = []
        for val in value:
            new_values.append(utils.round_digits(val))
        return ",".join(new_values)

class _DepthTag(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        super(self.__class__,
              self).__init__(MUTECT_ABBREVIATION,
                             common_tags.DEPTH_TAG,
                             'Jacquard depth for MuTect (based on DP)')

    def add_tag_values(self, vcf_record):
        if "DP" in vcf_record.format_tags:
            sample_values = {}
            for samp in vcf_record.sample_tag_values:
                sample_values[samp] = vcf_record.sample_tag_values[samp]["DP"]
            vcf_record.add_sample_tag_value(self.tag_id, sample_values)

class _SomaticTagSS(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods
    _DESCRIPTION = (\
'''Jacquard somatic status for MuTect: 0=non-somatic,1=somatic (based on SS
 FORMAT tag)''').replace("\n","")
    #pylint: disable=too-few-public-methods
    def __init__(self):
        super(self.__class__,
              self).__init__(MUTECT_ABBREVIATION,
                             common_tags.SOMATIC_TAG,
                             self._DESCRIPTION)

    def add_tag_values(self, vcf_record):
        sample_values = {}
        if "SS" in vcf_record.format_tags:
            for sample in vcf_record.sample_tag_values:
                somatic_status = vcf_record.sample_tag_values[sample]["SS"]
                sample_values[sample] = self._somatic_status(somatic_status)
        else:
            for sample in vcf_record.sample_tag_values:
                sample_values[sample] = "0"
        vcf_record.add_sample_tag_value(self.tag_id, sample_values)

    @staticmethod
    def _somatic_status(ss_value):
        if ss_value == "2":
            return "1"
        else:
            return "0"

class _SomaticTagFilterMutectCalls(common_tags.AbstractJacquardTag):
    #pylint: disable=too-few-public-methods
    _DESCRIPTION = (\
'''Jacquard somatic status for MuTect: 0=non-somatic,1=somatic (based on
 FilterMutectCalls setting filter to PASS)''').replace("\n","")
    #pylint: disable=too-few-public-methods
    def __init__(self):
        super(self.__class__,
              self).__init__(MUTECT_ABBREVIATION,
                             common_tags.SOMATIC_TAG,
                             self._DESCRIPTION)

    def add_tag_values(self, vcf_record):
        sample_values = {}
        if vcf_record.filter == "PASS":
            for sample in vcf_record.sample_tag_values:
                sample_values[sample] = self._somatic_status(vcf_record, sample)
        else:
            for sample in vcf_record.sample_tag_values:
                sample_values[sample] = "0"
        vcf_record.add_sample_tag_value(self.tag_id, sample_values)

    @staticmethod
    def _somatic_status(vcf_record, sample):
        tag_values = vcf_record.sample_tag_values[sample]
        try:
            gt = tag_values["GT"]
        except KeyError:
            msg_fmt = ('Cannot assign somatic status using FilterMutectCalls '
                       'when sample GT absent: '
                       '(CHROM:POS:REF:ALT={}:{}:{}:{})')
            msg = msg_fmt.format(vcf_record.chrom,
                                 vcf_record.pos,
                                 vcf_record.ref,
                                 vcf_record.alt)
            raise utils.JQException(msg)

        if gt == "0/0" or gt == '0|0':
            return "0"
        else:
            return "1"

def _build_somatic_tag(metaheaders):
    if "##source=FilterMutectCalls" in metaheaders:
        return _SomaticTagFilterMutectCalls()
    else:
        return _SomaticTagSS()

class _Mutect1Parser(object):
    _MUTECT1_METAHEADER_REGEX = re.compile('^##MuTect=')
    _MUTECT1_METAHEADER_DICT = re.compile('^##MuTect="(.*)"')
    _MUTECT1_DICT_KEY_NORMAL = 'normal_sample_name'
    _MUTECT1_DICT_KEY_TUMOR = 'tumor_sample_name'


    @staticmethod
    def is_mutect_metaheader(metaheaders):
        for metaheader in metaheaders:
            if _Mutect1Parser._MUTECT1_METAHEADER_REGEX.search(metaheader):
                return True
        return False

    @staticmethod
    def build_mutect_dict(metaheaders, normal_key, tumor_key):
        def get_mutect_header(metaheaders):
            for metaheader in metaheaders:
                match = _Mutect1Parser._MUTECT1_METAHEADER_DICT.search(metaheader)
                if match:
                    return match.group(1)
            return None

        mutect_dict = {}
        header = get_mutect_header(metaheaders)
        for token in header.split():
            try:
                key, value = token.split("=")
                if key == _Mutect1Parser._MUTECT1_DICT_KEY_NORMAL:
                    mutect_dict[normal_key] = value
                elif key == _Mutect1Parser._MUTECT1_DICT_KEY_TUMOR:
                    mutect_dict[tumor_key] = value
            except ValueError:
                pass
        return mutect_dict

class _Mutect2Parser(object):
    _MUTECT2_METAHEADER_REGEX = re.compile('^##GATKCommandLine.*?=<.*ID=Mu[tT]ect2')
    _MUTECT2_METAHEADER_COMMAND_REGEX = re.compile('^##GATKCommandLine.*?=<.*ID=Mu[tT]ect2.*CommandLine.*?="(.*?)"')
    _MUTECT2_METAHEADER_OLD_SAMPLE_REGEX = re.compile('^##SAMPLE=<ID=(.*?),SampleName=(.*?),.*')
    _MUTECT2_DICT_KEY_NORMAL = 'normal_sample'
    _MUTECT2_DICT_KEY_TUMOR = 'tumor_sample'
    _MUTECT2_METAHEADER_NEW_SAMPLE_REGEX = re.compile(
        '^##(' + _MUTECT2_DICT_KEY_NORMAL + '|' + _MUTECT2_DICT_KEY_TUMOR + ')=(.*)')


    @staticmethod
    def is_mutect_metaheader(metaheaders):
        for metaheader in metaheaders:
            if _Mutect2Parser._MUTECT2_METAHEADER_REGEX.search(metaheader):
                return True
        return False

    @staticmethod
    def _mutect_dict_from_new_sample_metalines(metaheaders, normal_key, tumor_key):
        mutect_keys = {_Mutect2Parser._MUTECT2_DICT_KEY_NORMAL: normal_key,
                       _Mutect2Parser._MUTECT2_DICT_KEY_TUMOR: tumor_key}
        mutect_dict = {}
        for metaheader in metaheaders:
            match = _Mutect2Parser._MUTECT2_METAHEADER_NEW_SAMPLE_REGEX.search(metaheader)
            if match:
                key = mutect_keys[match.group(1)]
                mutect_dict[key] = match.group(2)
        return mutect_dict

    @staticmethod
    def _mutect_dict_from_old_sample_metalines(metaheaders, normal_key, tumor_key):
        mutect_keys = {'NORMAL': normal_key, 'TUMOR': tumor_key}
        mutect_dict = {}
        for metaheader in metaheaders:
            match = _Mutect2Parser._MUTECT2_METAHEADER_OLD_SAMPLE_REGEX.search(metaheader)
            if match and match.group(1) in mutect_keys:
                key = mutect_keys[match.group(1)]
                mutect_dict[key] = match.group(2)
        return mutect_dict

    @staticmethod
    def _mutect_dict_from_command_line(metaheaders, normal_key, tumor_key):
        def get_mutect_header(metaheaders):
            for metaheader in metaheaders:
                match = _Mutect2Parser._MUTECT2_METAHEADER_COMMAND_REGEX.search(metaheader)
                if match:
                    return match.group(1)
            return None

        command_line_args = get_mutect_header(metaheaders).split()
        header_parser = argparse.ArgumentParser(add_help=False)
        header_parser.add_argument("--normal-sample")
        header_parser.add_argument("--tumor-sample")
        args, _ = header_parser.parse_known_args(command_line_args)
        mutect_dict = {}
        if args.normal_sample:
            mutect_dict[normal_key] = args.normal_sample
        if args.tumor_sample:
            mutect_dict[tumor_key] = args.tumor_sample
        return mutect_dict

    @staticmethod
    def build_mutect_dict(metaheaders, normal_key, tumor_key):
        normal_tumor_strategies = [
                _Mutect2Parser._mutect_dict_from_new_sample_metalines,
                _Mutect2Parser._mutect_dict_from_old_sample_metalines,
                _Mutect2Parser._mutect_dict_from_command_line,
                ]
        for strategy in normal_tumor_strategies:
            mutect_dict = strategy(metaheaders, normal_key, tumor_key)
            if mutect_dict:
                return mutect_dict
        #
        # mutect_dict = _Mutect2Parser._mutect_dict_from_new_sample_metalines(metaheaders, normal_key, tumor_key)
        # if not mutect_dict:
        #     mutect_dict = _Mutect2Parser._mutect_dict_from_old_sample_metalines(metaheaders, normal_key, tumor_key)
        # if not mutect_dict:
        #     mutect_dict = _Mutect2Parser._mutect_dict_from_command_line(metaheaders, normal_key, tumor_key)
        # return mutect_dict

def _get_mutect_parser(metaheaders):
    if _Mutect1Parser.is_mutect_metaheader(metaheaders):
        return _Mutect1Parser
    elif _Mutect2Parser.is_mutect_metaheader(metaheaders):
        return _Mutect2Parser
    return None

class Mutect(object):
    #pylint: disable=too-few-public-methods
    """Recognize and transform MuTect VCFs to standard Jacquard format.

    MuTect VCFs are blessedly compliant and straightforward to translate, with
    the following exception. The incoming header has the sample name values
    (derived from the input alignments). To play well with other callers like
    Strelka and VarScan, the sample headers are replaced with Normal and Tumor.
    """
    _NORMAL_SAMPLE_KEY = "normal_sample_name"
    _TUMOR_SAMPLE_KEY = "tumor_sample_name"

    def __init__(self):
        self.name = "MuTect"
        self.abbr = "MT"

    @staticmethod
    def _is_mutect_vcf(file_reader):
        if not file_reader.file_name.lower().endswith(".vcf"):
            return False
        vcf_reader = vcf.VcfReader(file_reader)
        return _get_mutect_parser(vcf_reader.metaheaders) != None

    def claim(self, file_readers):
        """Recognizes and claims MuTect VCFs form the set of all input VCFs.

        Each defined caller has a chance to evaluate and claim all the incoming
        files as something that it can process.

        Args:
            file_readers: the collection of currently unclaimed files

        Returns:
            A tuple of unclaimed readers and MuTectVcfReaders.
        """
        unclaimed_readers = []
        vcf_readers = []
        for file_reader in file_readers:
            if self._is_mutect_vcf(file_reader):
                vcf_reader = vcf.VcfReader(file_reader)
                vcf_readers.append(_MutectVcfReader(vcf_reader))
            else:
                unclaimed_readers.append(file_reader)
        return (unclaimed_readers, vcf_readers)

    @staticmethod
    def _build_mutect_dict(metaheaders):
        mutect_parser = _get_mutect_parser(metaheaders)
        return mutect_parser.build_mutect_dict(metaheaders,
                                               Mutect._NORMAL_SAMPLE_KEY,
                                               Mutect._TUMOR_SAMPLE_KEY)

    def _get_new_column_header(self, vcf_reader):
        """Returns a standardized column header.

        MuTect sample headers include the name of input alignment, which is
        nice, but doesn't match up with the sample names reported in Strelka
        or VarScan. To fix this, we replace with NORMAL and TUMOR using the
        MuTect metadata command line to replace them correctly."""
        mutect_dict = self._build_mutect_dict(vcf_reader.metaheaders)

        new_header_list = []
        required_keys = set([self._NORMAL_SAMPLE_KEY, self._TUMOR_SAMPLE_KEY])
        mutect_keys = set(mutect_dict.keys())

        if not required_keys.issubset(mutect_keys):
            raise utils.JQException("Unable to determine normal "
                                    "and tumor sample ordering "
                                    "based on MuTect metaheader.")

        for field_name in vcf_reader.column_header.split("\t"):
            if field_name == mutect_dict[self._NORMAL_SAMPLE_KEY]:
                field_name = "NORMAL"
            elif field_name == mutect_dict[self._TUMOR_SAMPLE_KEY]:
                field_name = "TUMOR"
            new_header_list.append(field_name)

        return "\t".join(new_header_list)


class _MutectVcfReader(object):
    """Adapter that presents a MuTect VCF as a VcfReader.

    This follows the VcfReader interface, delegating calls to the base
    VcfReader, adjusting metaheaders, column header, and individual
    variants as appropriate.

    See VcfReader for more info.
    """
    def __init__(self, vcf_reader):
        self._vcf_reader = vcf_reader
        self._caller = Mutect()
        self.tags = [common_tags.ReportedTag(MUTECT_ABBREVIATION),
                     common_tags.PassedTag(MUTECT_ABBREVIATION),
                     _AlleleFreqTag(vcf_reader.metaheaders),
                     _DepthTag(),
                     _build_somatic_tag(vcf_reader.metaheaders),
                     _GenotypeTag()]

    def _get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

    @property
    def caller_name(self):
        return self._caller.name

    @property
    def file_name(self):
        return self._vcf_reader.file_name

    def open(self):
        return self._vcf_reader.open()

    def close(self):
        return self._vcf_reader.close()

    @staticmethod
    def expected_file_format():
        return [""]

    @property
    def metaheaders(self):
        new_metaheaders = list(self._vcf_reader.metaheaders)
        new_metaheaders.extend(self._get_new_metaheaders())
        caller_metaheader = "##jacquard.translate.caller={0}".\
                format(self._caller.name)
        new_metaheaders.append(caller_metaheader)

        return new_metaheaders

    @property
    def column_header(self):
        return self._caller._get_new_column_header(self._vcf_reader)

    def tagged_vcf_records(self):
        for vcf_record in self._vcf_reader.vcf_records():
            yield vcf_record

    def vcf_records(self):
        for vcf_record in self._vcf_reader.vcf_records():
            yield self._add_tags(vcf_record)

    def _add_tags(self, vcf_record):
        for tag in self.tags:
            tag.add_tag_values(vcf_record)
        return vcf_record

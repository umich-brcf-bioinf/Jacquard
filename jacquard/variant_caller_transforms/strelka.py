"""Interprets Strelka VCF files adding Jacquard standard information.

* Strelka VCFs are assumed to have a ".vcf" extension and have a
    "##source=strelka" metaheader.
* Strelka produces a separate file for SNVs and indels. Jacquard can process
    either or both.
* Jacquard standard tags are based on tier 2 data and use different source tags
    based on whether the file is indel or snp.

See tag definitions for more info.
"""
from __future__ import print_function, absolute_import, division

from collections import defaultdict, OrderedDict
import os
import re

import jacquard.utils.logger as logger
import jacquard.utils.utils as utils
import jacquard.variant_caller_transforms.common_tags as common_tags
import jacquard.utils.vcf as vcf


JQ_STRELKA_TAG = "JQ_SK_"
VERSION = "v2.0.15"

class _AlleleFreqTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}AF,'
                           'Number=A,'
                           'Type=Float,'
                           'Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth. Uses (TIR tier 2)/DP2 if available, otherwise uses (ACGT tier2 depth) / DP2)">')\
                           .format(JQ_STRELKA_TAG)

    @staticmethod
    def _get_tier2_base_depth(sample_format_dict, alt_allele):
        numerator = float(sample_format_dict[alt_allele + "U"].split(",")[1])
        tags = ["AU", "CU", "TU", "GU"]
        depth = 0

        for tag in tags:
            depth += float(sample_format_dict[tag].split(",")[1])
        freq = numerator/depth if depth != 0 else 0.0

        return freq

    def _get_snp_allele_freq_per_sample(self, vcf_record, sample):
        afs = []
        split_alt = vcf_record.alt.split(",")

        for alt_allele in split_alt:
            sample_format_dict = vcf_record.sample_tag_values[sample]
            freq = self._get_tier2_base_depth(sample_format_dict, alt_allele)

            rounded_af = self._standardize_af(str(freq))
            capped_af = min(rounded_af, "1.00")
            afs.append(capped_af)

        return afs

    def _get_indelallelefreq_per_sample(self, vcf_record, sample):
        afs = []
        indel_depths = vcf_record.sample_tag_values[sample]["TIR"]
        numerator = float(indel_depths.split(",")[1])
        denominator = float(vcf_record.sample_tag_values[sample]["DP2"])
        freq = numerator/denominator if denominator != 0 else 0.0

        rounded_af = self._standardize_af(str(freq))
        capped_af = min(rounded_af, "1.00")
        afs.append(capped_af)

        return afs

    def add_tag_values(self, vcf_record):
        sample_values = {}
        if vcf_record.alt == ".":
            for sample in vcf_record.sample_tag_values:
                sample_values[sample] = "."
        else:
            for sample in vcf_record.sample_tag_values:
                if "AU" in vcf_record.format_tags:#if it's an snp
                    afs = self._get_snp_allele_freq_per_sample(vcf_record,
                                                               sample)
                    sample_values[sample] = ",".join(afs)
                elif "TIR" in vcf_record.format_tags: #if it's an indel
                    afs = self._get_indelallelefreq_per_sample(vcf_record,
                                                               sample)
                    sample_values[sample] = ",".join(afs)
                else:
                    continue

        if sample_values:
            vcf_record.add_sample_tag_value(JQ_STRELKA_TAG + "AF",
                                            sample_values)

    @staticmethod
    def _standardize_af(value):
        return utils.round_two_digits(value)

class _DepthTag(object):
    #pylint: disable=too-few-public-methods
    REQUIRED_TAGS = set(["DP2", "AU"])
    NUCLEOTIDE_DEPTH_TAGS = ["AU", "CU", "TU", "GU"]

    @classmethod
    def _get_tier2_base_depth(cls, sample_tags):
        if "DP2" in sample_tags:
            return sample_tags["DP2"]
        depth = 0
        for tag in _DepthTag.NUCLEOTIDE_DEPTH_TAGS:
            tier2_depth = sample_tags[tag].split(",")[1]
            depth += int(tier2_depth)
        return str(depth)

    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}DP,'
                           'Number=1,'
                           'Type=Float,'
                           'Description="Jacquard depth for Strelka (uses DP2 if available, otherwise uses ACGT tier2 depth)">')\
                           .format(JQ_STRELKA_TAG)

    @staticmethod
    def add_tag_values(vcf_record):
        if vcf_record.format_tags.isdisjoint(_DepthTag.REQUIRED_TAGS):
            return
        sample_values = {}
        for sample in vcf_record.sample_tag_values:
            sample_tags = vcf_record.sample_tag_values[sample]
            sample_values[sample] = _DepthTag._get_tier2_base_depth(sample_tags)
        vcf_record.add_sample_tag_value(JQ_STRELKA_TAG + "DP", sample_values)


##TODO (cgates): Make this robust to sample order changes
class _SomaticTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}HC_SOM,'
                           'Number=1,Type=Integer,'
                           'Description="Jacquard somatic status for Strelka: 0=non-somatic,1=somatic (based on PASS in FILTER column)">')\
                           .format(JQ_STRELKA_TAG)

    @staticmethod
    def add_tag_values(vcf_record):
        sample_values = {}
        for i, sample in enumerate(vcf_record.sample_tag_values):
            sample_values[sample] = _SomaticTag._somatic_status(i, vcf_record)
        strelka_tag = JQ_STRELKA_TAG + "HC_SOM"
        vcf_record.add_sample_tag_value(strelka_tag, sample_values)

    @classmethod
    def _somatic_status(cls, sample_index, vcf_record):
        if sample_index == 1 and vcf_record.filter == "PASS":
            return "1"
        else:
            return "0"

class Strelka(object):
    """Recognize and transform Strelka VCFs to standard Jacquard format.

    Note that Strelka sometimes reports variants which fail the filter as
    having an ALT of "."; this seems to be Strelka's shorthand for saying,
    "I couldn't be sure there was an ALT". Unfortunately, that's not a valid
    ALT value so these rows are flagged and (in a typical workflow) excluded.
    """

    def __init__(self):
        self.name = "Strelka"
        self.abbr = "SK"
        self.meta_header = "##jacquard.normalize_strelka.sources={0},{1}\n"

    ##TODO (cgates): deprecated; remove
    @staticmethod
    def validate_input_file(meta_headers, dummy_column_header):
        return "##source=strelka" in meta_headers

    @staticmethod
    def _is_strelka_vcf(file_reader):
        if file_reader.file_name.endswith(".vcf"):
            vcf_reader = vcf.VcfReader(file_reader)
            return "##source=strelka" in vcf_reader.metaheaders
        return False

    def _find_strelka_files(self, file_readers):
        prefix_dict = OrderedDict()
        unclaimed_readers = []
        for file_reader in file_readers:
            if self._is_strelka_vcf(file_reader):
                prefix, _ = os.path.splitext(file_reader.file_name)
                prefix_dict[prefix] = file_reader
            else:
                unclaimed_readers.append(file_reader)

        return prefix_dict, unclaimed_readers

    @staticmethod
    def _split_prefix_by_patient(prefix_dict):
        prefix_by_patients = defaultdict(list)
        snv_indel = ["snvs", "indels"]

        for prefix in prefix_dict:
            patient_names = [i for i in prefix.split(".") if i not in snv_indel]
            patient_name = ".".join(patient_names)
            prefix_by_patients[patient_name].append(prefix)

        return prefix_by_patients

    @staticmethod
    def _validate_vcf_readers(prefix_by_patients):
        number_of_files = set()
        for file_names in prefix_by_patients.values():
            if len(file_names) == 1:
                for file_name in file_names:
                    if re.search("snvs", file_name):
                        msg = "Strelka VCF [{}] has no indels file."
                        logger.error(msg, file_name)
                    elif re.search("indels", file_name):
                        msg = "Strelka VCF [{}] has no snvs file."
                        logger.error(msg, file_name)
            number_of_files.add(len(file_names))

        if len(number_of_files) > 1:
            msg = ("Some Strelka VCFs were missing either a snvs or indels "
                   "file. Review inputs/command options and try again.")
            raise utils.JQException(msg)

    @staticmethod
    def _create_vcf_readers(prefix_to_readers):
        vcf_readers = []
        for file_reader in prefix_to_readers.values():
            vcf_reader = vcf.VcfReader(file_reader)
            vcf_readers.append(_StrelkaVcfReader(vcf_reader))

        return vcf_readers

    def claim(self, file_readers):
        """Recognizes and claims Strelka VCFs form the set of all input VCFs.

        Each defined caller has a chance to evaluate and claim all the incoming
        files as something that it can process.

        Args:
            file_readers: the collection of currently unclaimed files

        Returns:
            A tuple of unclaimed readers and StrelkaVcfReaders.
        """
        (prefix_to_reader,
         unclaimed_readers) = self._find_strelka_files(file_readers)
        prefix_by_patients = self._split_prefix_by_patient(prefix_to_reader)
        self._validate_vcf_readers(prefix_by_patients)
        vcf_readers = self._create_vcf_readers(prefix_to_reader)

        return (unclaimed_readers, vcf_readers)

class _StrelkaVcfReader(object):
    """Adapter that presents a Strelka VCF as a VcfReader.

    This follows the VcfReader interface, delegating calls to the base
    VcfReader, adjusting metaheaders, and individual
    variants as appropriate.

    See VcfReader for more info.
    """
    def __init__(self, vcf_reader):
        self._vcf_reader = vcf_reader
        self._caller = Strelka()
        self.tags = [common_tags.ReportedTag(JQ_STRELKA_TAG),
                     common_tags.PassedTag(JQ_STRELKA_TAG),
                     _AlleleFreqTag(),
                     _DepthTag(),
                     _SomaticTag()]

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
        return ["snvs", "indels"]

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
        return self._vcf_reader.column_header

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


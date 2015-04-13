"""Interprets MuTect VCF files adding Jacquard standard information.

MuTect VCFs are assumed to have a ".vcf" extension and have a valid
"##MuTect=..." metaheader.
"""
from __future__ import print_function, absolute_import, division
import jacquard.utils.utils as utils
import jacquard.variant_caller_transforms.common_tags as common_tags
import jacquard.utils.vcf as vcf

JQ_MUTECT_TAG = "JQ_MT_"
VERSION = "v1.1.4"

class _AlleleFreqTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = ('##FORMAT=<ID={0}AF,'
                           'Number=A,'
                           'Type=Float,'
                           #pylint: disable=line-too-long
                           'Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA)">')\
                           .format(JQ_MUTECT_TAG)

    def add_tag_values(self, vcf_record):
        if "FA" in vcf_record.format_tags:
            sample_values = {}
            for sample in vcf_record.sample_tag_values:
                freq = vcf_record.sample_tag_values[sample]["FA"].split(",")
                sample_values[sample] = self._standardize_af(freq)
            vcf_record.add_sample_tag_value(JQ_MUTECT_TAG + "AF", sample_values)

    @staticmethod
    def _standardize_af(value):
        new_values = []
        for val in value:
            new_values.append(utils.round_two_digits(val))
        return ",".join(new_values)

class _DepthTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = ('##FORMAT=<ID={0}DP,'
                           'Number=1,'
                           'Type=Float,'
                           #pylint: disable=line-too-long
                           'Description="Jacquard depth for MuTect (based on DP)">')\
                           .format(JQ_MUTECT_TAG)

    @staticmethod
    def add_tag_values(vcf_record):
        if "DP" in vcf_record.format_tags:
            sample_values = {}
            for samp in vcf_record.sample_tag_values:
                sample_values[samp] = vcf_record.sample_tag_values[samp]["DP"]
            vcf_record.add_sample_tag_value(JQ_MUTECT_TAG + "DP", sample_values)

class _SomaticTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.metaheader = ('##FORMAT=<ID={0}HC_SOM,'
                           'Number=1,'
                           'Type=Integer,'
                           #pylint: disable=line-too-long
                           'Description="Jacquard somatic status for MuTect: 0=non-somatic,1=somatic (based on SS FORMAT tag)">')\
                           .format(JQ_MUTECT_TAG)

    def add_tag_values(self, vcf_record):
        mutect_tag = JQ_MUTECT_TAG + "HC_SOM"
        sample_values = {}
        if "SS" in vcf_record.format_tags:
            for sample in vcf_record.sample_tag_values:
                somatic_status = vcf_record.sample_tag_values[sample]["SS"]
                sample_values[sample] = self._somatic_status(somatic_status)
        else:
            for sample in vcf_record.sample_tag_values:
                sample_values[sample] = "0"
        vcf_record.add_sample_tag_value(mutect_tag, sample_values)

    @staticmethod
    def _somatic_status(ss_value):
        if ss_value == "2":
            return "1"
        else:
            return "0"

class Mutect(object):
    """Recognize and transform MuTect VCFs to standard Jacquard format.

    MuTect VCFs are blessedly compliant and straightforward to translate, with
    the following exception. The incoming header has the sample name values
    (derived from the input alignments). To play well with other callers like
    Strelka and VarScan, the sample headers are replaced with Normal and Tumor.
    """
    _MUTECT_METAHEADER_PREFIX = "##MuTect"
    _NORMAL_SAMPLE_KEY = "normal_sample_name"
    _TUMOR_SAMPLE_KEY = "tumor_sample_name"

    def __init__(self):
        self.name = "MuTect"
        self.abbr = "MT"

    ##TODO (cgates): deprecate; remove
    @staticmethod
    def validate_input_file(meta_headers, dummy_column_header):
        for line in meta_headers:
            if line.startswith(Mutect._MUTECT_METAHEADER_PREFIX):
                return True
        return False

    @staticmethod
    def _is_mutect_vcf(file_reader):
        if file_reader.file_name.lower().endswith(".vcf"):
            vcf_reader = vcf.VcfReader(file_reader)
            for metaheader in vcf_reader.metaheaders:
                if metaheader.startswith(Mutect._MUTECT_METAHEADER_PREFIX):
                    return True
        return False

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
        mutect_dict = {}
        for metaheader in metaheaders:
            if metaheader.startswith(Mutect._MUTECT_METAHEADER_PREFIX):
                split_line = metaheader.strip('"').split(" ")
                for item in split_line:
                    split_item = item.split("=")
                    try:
                        mutect_dict[split_item[0]] = split_item[1]
                    except IndexError:
                        pass

        return mutect_dict

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
        self.tags = [common_tags.ReportedTag(JQ_MUTECT_TAG),
                     common_tags.PassedTag(JQ_MUTECT_TAG),
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


"""Interprets VarScan2 VCF files adding Jacquard standard information.

* VarScan VCFs are assumed to have a ".vcf" extension and have a
    "##source=VarScan2" metaheader.
* VarScan produces a separate file for SNPs and indels. Jacquard can process
    either or both.
* Jacquard requires the VarScan VCF outputs. The VarScan workflow has optional
    extra steps to:
    a) partition the results into Germline, LOH, and Somatic files
    b) filter to a subset of high-confidence variants
    If these files are provided in the input, Jacquard will use them to
    flag variants not found in the high-confidence files so that they can be
    optionally filtered out.
* If provided, high-confidence files should follow the naming convention of
    patientName.*.fpfilter.pass
    The high-confidence file suffix can be supplied as a command line arg.

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


_VARSCAN_SOMATIC_HEADER = ("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|"
                           "NORMAL|TUMOR").replace("|", "\t")
JQ_VARSCAN_TAG = "JQ_VS_"
VERSION = "v2.3"

class _AlleleFreqTag(object):
    #pylint: disable=too-few-public-methods
    @classmethod
    def _standardize_af(cls, value):
        new_values = []
        for val in value:
            new_val = str(float(val.strip("%")) / 100)
            new_values.append(utils.round_two_digits(new_val))
        return ",".join(new_values)

    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}AF,'
                           'Number=A,'
                           'Type=Float,'
                           'Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ)">')\
                           .format(JQ_VARSCAN_TAG)

    @staticmethod
    def add_tag_values(vcf_record):
        sample_values = {}
        if "FREQ" in vcf_record.format_tags:
            for sample in vcf_record.sample_tag_values:
                freq = vcf_record.sample_tag_values[sample]["FREQ"].split(",")
                sample_values[sample] = _AlleleFreqTag._standardize_af(freq)
            vcf_record.add_sample_tag_value(JQ_VARSCAN_TAG + "AF",
                                            sample_values)

class _DepthTag(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}DP,'
                           'Number=1,'
                           'Type=Float,'
                           'Description="Jacquard depth for VarScan (based on DP)">')\
                           .format(JQ_VARSCAN_TAG)

    @staticmethod
    def add_tag_values(vcf_record):
        if "DP" in vcf_record.format_tags:
            sample_values = {}
            for sample in vcf_record.sample_tag_values:
                depth = vcf_record.sample_tag_values[sample]["DP"]
                sample_values[sample] = depth
            vcf_record.add_sample_tag_value(JQ_VARSCAN_TAG + "DP",
                                            sample_values)

##TODO (cgates): Make this robust to sample order changes
class _SomaticTag(object):
    #pylint: disable=too-few-public-methods

    @staticmethod
    def _somatic_status(sample_index):
        if sample_index == 0:  # it's NORMAL
            return "0"
        else:  # it's TUMOR
            return "1"

    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}HC_SOM,'
                           'Number=1,'
                           'Type=Integer,'
                           'Description="Jacquard somatic status for VarScan: 0=non-somatic,1=somatic (based on SOMATIC info tag and if sample is TUMOR)">')\
                           .format(JQ_VARSCAN_TAG)

    @staticmethod
    def add_tag_values(vcf_record):
        info_array = vcf_record.info.split(";")
        varscan_tag = JQ_VARSCAN_TAG + "HC_SOM"
        sample_values = {}

        if "SS=2" in info_array and vcf_record.filter == "PASS":
            for i, sample in enumerate(vcf_record.sample_tag_values):
                sample_values[sample] = _SomaticTag._somatic_status(i)
        else:
            for sample in vcf_record.sample_tag_values:
                sample_values[sample] = "0"

        vcf_record.add_sample_tag_value(varscan_tag, sample_values)

#TODO: (cgates/jebene): All tags should have _TAG_ID as implemented below
class _HCTag(object):
    #pylint: disable=too-few-public-methods
    _FILTERS_TO_REPLACE = set(["", ".", "pass"])
    _TAG_ID = "{}LOW_CONFIDENCE".format(JQ_VARSCAN_TAG)

    def __init__(self, file_reader):
        #pylint: disable=line-too-long
        self.metaheader = ('##FILTER=<ID={},'
                           'Number=1,'
                           'Type=Flag,'
                           'Description="Jacquard high-confidence somatic flag for VarScan. Based on intersection with filtered VarScan variants">')\
                           .format(self._TAG_ID)
        self._hc_loci = self._parse_file_reader(file_reader)

    @staticmethod
    def _parse_file_reader(file_reader):
        column_header = None
        hc_loci = set()
        file_reader.open()
        for line in file_reader.read_lines():
            if line.startswith("chrom\tposition"):
                column_header = line
            else:
                split_line = line.split("\t")
                hc_loci.add((split_line[0], split_line[1]))
        file_reader.close()

        #TODO : (cgates): Please test this
        if not column_header:
            raise utils.JQException("Error. The hc file {} is in an incorrect"
                                    "format. Review inputs and try"
                                    "again.".format(file_reader.file_name))

        return hc_loci

    def add_tag_values(self, vcf_record):
        if (vcf_record.chrom, vcf_record.pos) not in self._hc_loci:
            vcf_record.add_or_replace_filter(self._TAG_ID)
        return vcf_record

class Varscan(object):
    """Recognize and transform VarScan VCFs to standard Jacquard format."""

    _DEFAULT_REGEX = "Somatic.hc.fpfilter.pass"

    def __init__(self, args=None):
        self.name = "VarScan"
        self.abbr = "VS"
        self.meta_header = "##jacquard.normalize_varscan.sources={0},{1}\n"

        self.hc_file_pattern = self._get_hc_file_pattern(args)

    @staticmethod
    def _get_hc_file_pattern(args):
        msg = ("The specified regex [{}] could not be compiled. "
               "Review inputs and try again")
        if args and args.varscan_hc_filter_file_regex:
            try:
                compiled_regex = re.compile(args.varscan_hc_filter_file_regex)
            except:
                raise utils.UsageError(msg, args.varscan_hc_filter_file_regex)
        else:
            try:
                compiled_regex = re.compile(Varscan._DEFAULT_REGEX)
            except:
                raise utils.UsageError(msg, Varscan._DEFAULT_REGEX)

        return compiled_regex

    ##TODO (cgates): deprecated; remove
    @staticmethod
    def validate_input_file(meta_headers, column_header):
        if "##source=VarScan2" not in meta_headers:
            return 0

        if _VARSCAN_SOMATIC_HEADER == column_header:
            return 1
        else:
            raise utils.JQException("Unexpected VarScan VCF structure - "
                                    "missing NORMAL and TUMOR headers.")

    @staticmethod
    def _validate_filter_file(file_reader):
        column_header = 0
        file_reader.open()
        for line in file_reader.read_lines():
            if line.startswith("chrom\tposition"):
                column_header = line
                break
        file_reader.close()

        if column_header:
            return file_reader

    @staticmethod
    def _is_varscan_vcf(file_reader):
        if file_reader.file_name.endswith(".vcf"):
            vcf_reader = vcf.VcfReader(file_reader)
            return "##source=VarScan2" in vcf_reader.metaheaders
        return False

#TODO: (cgates) Add check of header line - extract constant from HCTag?
    def _is_varscan_hc_file(self, file_reader):
        return self.hc_file_pattern.search(file_reader.file_name)

    @staticmethod
    def _raise_invalid_filter_exception(invalid_filter_files):
        if invalid_filter_files[5:]:
            omitted_files = ("...({} file(s) omitted)")\
                            .format(len(invalid_filter_files[5:]))
        else:
            omitted_files = ""
        first_five_fnames = [i.file_name for i in invalid_filter_files[:5]]
        raise utils.JQException("The [{}] input files [{}{}] match "
                                "high-confidence file names, but the "
                                "file header is invalid or missing. "
                                "Review inputs and try again.",
                                len(invalid_filter_files),
                                first_five_fnames,
                                omitted_files)

    @staticmethod
    def _validate_file_pairing(unpaired_hc_files, unpaired_vcf_files):
        if unpaired_hc_files:
            for unpaired_hc_file in unpaired_hc_files[:5]:
                msg = ("The VarScan high-confidence file [{}] has no matching "
                       "VCF file.")
                logger.error(msg, unpaired_hc_file.file_name)

            msg = ("[{}] VarScan high-confidence file(s) did not have a "
                   "matching VCF file. See log for more details. Review inputs "
                   "and try again.")
            raise utils.JQException(msg, len(unpaired_hc_files))

        if unpaired_vcf_files:
            for unpaired_vcf_file in unpaired_vcf_files[:5]:
                msg = ("The VarScan VCF file [{}] has no matching "
                       "high-confidence file.")
                logger.error(msg, unpaired_vcf_file.file_name)

            msg = ("[{}] VarScan VCF file(s) did not have a matching "
                   "high-confidence file. See log for more details. Review "
                   "inputs and try again.")
            raise utils.JQException(msg, len(unpaired_vcf_files))

    def _find_varscan_files(self, file_readers):
        unclaimed_set = set()
        prefix_file_readers = OrderedDict()
        filter_files = set()

        for file_reader in file_readers:
            if self._is_varscan_vcf(file_reader):
                prefix, _ = os.path.splitext(file_reader.file_name)
                prefix_file_readers[prefix] = file_reader
            elif self._is_varscan_hc_file(file_reader):
                filter_files.add(file_reader)
            else:
                unclaimed_set.add(file_reader)

        if not filter_files:
            if not self.hc_file_pattern.match(Varscan._DEFAULT_REGEX):
                msg = ("The VarScan high-confidence filename regex [{}] "
                       "didn't match any files in the input directory. "
                       "The beginning of the high-confidence filename must "
                       "exactly match a VCF filename up to the .vcf extention. "
                       "Review inputs/command options and try again.")
                raise utils.UsageError(msg, self.hc_file_pattern.pattern)

        return prefix_file_readers, filter_files, unclaimed_set

    @staticmethod
    def _split_prefix_by_patient(prefix_vcf_readers):
        prefix_by_patients = defaultdict(list)
        snp_indel = ["snp", "indel"]
        for prefix in prefix_vcf_readers:
            patient_names = [i for i in prefix.split(".") if i not in snp_indel]
            patient_name = ".".join(patient_names)
            prefix_by_patients[patient_name].append(prefix)

        return prefix_by_patients

    @staticmethod
    def _validate_vcf_readers(prefix_by_patients):
        number_of_files = set()
        for file_names in prefix_by_patients.values():
            if len(file_names) == 1:
                for file_name in file_names:
                    if re.search("snp", file_name):
                        msg = "VarScan VCF [{}] has no indel file."
                        logger.error(msg, file_name)
                    elif re.search("indel", file_name):
                        msg = "VarScan VCF [{}] has no snp file."
                        logger.error(msg, file_name)
            number_of_files.add(len(file_names))

        if len(number_of_files) > 1:
            msg = ("Some Varscan VCFs were missing either a snp or indel file. "
                   "Review inputs/command options and try again.")
            raise utils.JQException(msg)

    @staticmethod
    def _remove_paired_filters(vcf_dict, filter_dict):
        for filter_file in vcf_dict.values():
            try:
                filter_dict.pop(filter_file)
            except KeyError:
                pass
        return filter_dict

    @staticmethod
    def _dictionaries_to_tuples(vcf_dict, filter_dict):
        flipped_filter = OrderedDict((y, x) for x, y in filter_dict.items())
        pairs = vcf_dict.copy()
        pairs.update(flipped_filter)
        return [(k, v) for k, v in pairs.items()]

    def _pair_files(self, prefixes, filter_files):
        vcf_dict = OrderedDict()
        filter_dict = OrderedDict()
        invalid_filter_files = []

        for prefix, file_reader in prefixes.items():
            vcf_dict[file_reader] = None
            for filter_reader in filter_files:
                filter_dict[filter_reader] = None

                if filter_reader.file_name.startswith(prefix):
                    if self._validate_filter_file(filter_reader):
                        vcf_dict[file_reader] = filter_reader
                    else:
                        invalid_filter_files.append(filter_reader)

        if invalid_filter_files:
            self._raise_invalid_filter_exception(invalid_filter_files)
        filter_dict = self._remove_paired_filters(vcf_dict, filter_dict)

        return self._dictionaries_to_tuples(vcf_dict, filter_dict)

    def _validate_file_pairs(self, tuples):
        unpaired_vcf_files = []
        unpaired_hc_files = []

        hc_values = set([i[1] for i in tuples])
        if len(hc_values) != 1:
            for vcf_file_reader, hc_file_reader in tuples:
                if vcf_file_reader and not hc_file_reader:
                    unpaired_vcf_files.append(vcf_file_reader)
                if hc_file_reader and not vcf_file_reader:
                    unpaired_hc_files.append(hc_file_reader)

        self._validate_file_pairing(unpaired_hc_files, unpaired_vcf_files)

    @staticmethod
    def _create_vcf_readers(pair_tuples):
        vcf_readers = []
        for vcf_file_reader, hc_file_reader in pair_tuples:
            if vcf_file_reader and hc_file_reader:
                vcf_reader = _VarscanVcfReader(vcf.VcfReader(vcf_file_reader),
                                               hc_file_reader)
            elif vcf_file_reader and not hc_file_reader:
                vcf_reader = _VarscanVcfReader(vcf.VcfReader(vcf_file_reader))
            vcf_readers.append(vcf_reader)

        return vcf_readers

#pylint: disable=too-many-locals
    def claim(self, file_readers):
        """Recognizes and claims MuTect VCFs form the set of all input VCFs.

        Each defined caller has a chance to evaluate and claim all the incoming
        files as something that it can process. Since VarScan can claim
        high-confidence files as well, this process is significantly more
        complex than for other callers.

        Args:
            file_readers: the collection of currently unclaimed files

        Returns:
            A tuple of unclaimed readers and MuTectVcfReaders.
        """

        (prefix_to_readers,
         filter_files,
         unclaimed_set) = self._find_varscan_files(file_readers)

        prefix_by_patients = self._split_prefix_by_patient(prefix_to_readers)
        self._validate_vcf_readers(prefix_by_patients)
        vcf_filter_tuples = self._pair_files(prefix_to_readers, filter_files)
        self._validate_file_pairs(vcf_filter_tuples)
        vcf_readers = self._create_vcf_readers(vcf_filter_tuples)

        return list(unclaimed_set), vcf_readers

class _VarscanVcfReader(object):
    """Adapter that presents a VarScan VCF as a VcfReader.

    This follows the VcfReader interface, delegating calls to the base
    VcfReader, adjusting metaheaders and individual
    variants as appropriate.

    See VcfReader for more info.
    """
    def __init__(self, vcf_reader, som_hc_file_reader=None):
        self._vcf_reader = vcf_reader
        self._som_hc_file_reader = som_hc_file_reader
        self._caller = Varscan()
        self.tags = [common_tags.ReportedTag(JQ_VARSCAN_TAG),
                     common_tags.PassedTag(JQ_VARSCAN_TAG),
                     _AlleleFreqTag(),
                     _DepthTag(),
                     _SomaticTag()]

        if som_hc_file_reader:
            self.tags.insert(0, _HCTag(som_hc_file_reader))

    @property
    def file_name(self):
        return self._vcf_reader.file_name

    def _get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

    @property
    def caller_name(self):
        return self._caller.name

    def open(self):
        return self._vcf_reader.open()

    def close(self):
        return self._vcf_reader.close()

    @staticmethod
    def expected_file_format():
        return ["snp", "indel"]

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

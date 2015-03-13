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
from __future__ import print_function, absolute_import

from collections import defaultdict, OrderedDict
import os
import re

import jacquard.utils as utils
import jacquard.variant_callers.common_tags as common_tags
import jacquard.vcf as vcf


_VARSCAN_SOMATIC_HEADER = ("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|"
                           "NORMAL|TUMOR").replace("|", "\t")
JQ_VARSCAN_TAG = "JQ_VS_"


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

    def __init__(self):
        self.name = "VarScan"
        self.abbr = "VS"
        self.meta_header = "##jacquard.normalize_varscan.sources={0},{1}\n"
        self.hc_file_pattern = re.compile("fpfilter.pass")

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
    def _is_varscan_vcf(file_reader):
        if file_reader.file_name.endswith(".vcf"):
            vcf_reader = vcf.VcfReader(file_reader)
            return "##source=VarScan2" in vcf_reader.metaheaders
        return False

    #TODO: (cgates): Add check of header line (extract constant from HCTag?)
    def _is_varscan_hc_file(self, file_reader):
        return self.hc_file_pattern.search(file_reader.file_name)

    @staticmethod
    def _get_files_per_patient(file_readers):
        patient_to_files = defaultdict(list)
        for file_reader in sorted(file_readers):
            filename = file_reader.file_name
            patient = filename.split(".")[0]
            patient_to_files[patient].append(file_reader)

        return patient_to_files

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
        files_per_patient = self._get_files_per_patient(file_readers)

        unclaimed_set = set()
        trans_vcf_readers = []

        for patient in files_per_patient:
            prefix_reader = OrderedDict()
            filter_files = set()
            for file_reader in files_per_patient[patient]:
                if self._is_varscan_vcf(file_reader):
                    prefix, _ = os.path.splitext(file_reader.file_name)
                    prefix_reader[prefix] = file_reader
                elif self._is_varscan_hc_file(file_reader):
                    filter_files.add(file_reader)
                else:
                    unclaimed_set.add(file_reader)

            for prefix, reader in prefix_reader.items():
                vcf_reader = _VarscanVcfReader(vcf.VcfReader(reader))
                for filter_file in list(filter_files):
                    if filter_file.file_name.startswith(prefix):
                        vcf_reader = _VarscanVcfReader(vcf.VcfReader(reader),
                                                       filter_file)
                        filter_files.remove(filter_file)

                trans_vcf_readers.append(vcf_reader)
            unclaimed_set.update(filter_files)

        return list(unclaimed_set), trans_vcf_readers


#TODO: (cgates): If we can, I would rather inflate the high confidence set when
# we open and not on construction. There is a pretty safe/clean way to do this.
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

    def vcf_records(self):
        for vcf_record in self._vcf_reader.vcf_records():
            yield self._add_tags(vcf_record)

    def _add_tags(self, vcf_record):
        for tag in self.tags:
            tag.add_tag_values(vcf_record)
        return vcf_record

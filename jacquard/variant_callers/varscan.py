from __future__ import print_function, absolute_import
import jacquard.vcf as vcf
import jacquard.variant_callers.common_tags as common_tags
import jacquard.utils as utils
from jacquard import __version__
import os
import re


VARSCAN_SOMATIC_HEADER = ("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|"
                          "NORMAL|TUMOR").replace("|", "\t")
JQ_VARSCAN_TAG = "JQ_VS_"
_LOW_CONFIDENCE_FILTER = "JQ_LOW_CONFIDENCE"

class _AlleleFreqTag(object):
    #pylint: disable=too-few-public-methods
    @classmethod
    def _standardize_af(cls, value):
        new_values = []
        for val in value:
            new_val = str(float(val.strip("%"))/100)
            new_values.append(utils.round_two_digits(new_val))
        return ",".join(new_values)

    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}AF,'
                           'Number=A,'
                           'Type=Float,'
                           'Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ)",'
                           'Source="Jacquard",'
                           'Version={1}>').format(JQ_VARSCAN_TAG,
                                                  __version__)

    @staticmethod
    def add_tag_values(vcf_record):
        sample_values = {}
        if "FREQ" in vcf_record.format_tags:
            for sample in vcf_record.sample_tag_values:
                freq = vcf_record.sample_tag_values[sample]["FREQ"].split(",")
                sample_values[sample] = _AlleleFreqTag._standardize_af(freq)
            vcf_record.add_sample_tag_value(JQ_VARSCAN_TAG + "AF",
                                            sample_values)

#pylint: disable=too-few-public-methods
class _DepthTag(object):
    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}DP,'
                           'Number=1,'
                           'Type=Float,'
                           'Description="Jacquard depth for VarScan (based on DP)",'
                           'Source="Jacquard",'
                           'Version={1}>').format(JQ_VARSCAN_TAG,
                                                  __version__)

    @staticmethod
    def add_tag_values(vcf_record):
        if "DP" in vcf_record.format_tags:
            sample_values = {}
            for sample in vcf_record.sample_tag_values:
                depth = vcf_record.sample_tag_values[sample]["DP"]
                sample_values[sample] = depth
            vcf_record.add_sample_tag_value(JQ_VARSCAN_TAG + "DP",
                                            sample_values)

#pylint: disable=too-few-public-methods
class _SomaticTag(object):
    @staticmethod
    def _somatic_status(sample_index):
        if sample_index == 0: #it's NORMAL
            return "0"
        else: #it's TUMOR
            return "1"

    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}HC_SOM,'
                           'Number=1,'
                           'Type=Integer,'
                           'Description="Jacquard somatic status for VarScan: 0=non-somatic,1=somatic (based on SOMATIC info tag and if sample is TUMOR)",'
                           'Source="Jacquard",'
                           'Version={1}>').format(JQ_VARSCAN_TAG,
                                                  __version__)

    @staticmethod
    #TODO: jebene - edit this to reflect new HCTag changes
    def add_tag_values(vcf_record):
        info_array = vcf_record.info.split(";")
        varscan_tag = JQ_VARSCAN_TAG + "HC_SOM"
        sample_values = {}
        if "SS=2" in info_array and JQ_VARSCAN_TAG + "HC" in info_array:
            for i, sample in enumerate(vcf_record.sample_tag_values):
                sample_values[sample] = _SomaticTag._somatic_status(i)
        else:
            for sample in vcf_record.sample_tag_values:
                sample_values[sample] = "0"

        vcf_record.add_sample_tag_value(varscan_tag, sample_values)

class _HCTag(object):
    def __init__(self, file_reader):
        #pylint: disable=line-too-long
        self.metaheader = ('##FILTER=<ID={}HC,'
                           'Number=1,'
                           'Type=Flag,'
                           'Description="Jacquard high-confidence somatic flag for VarScan. Based on intersection with filtered VarScan variants"'
                           'Source="Jacquard",'
                           'Version={}>').format(JQ_VARSCAN_TAG,
                                                 __version__)
        self.hc_loci = self._parse_file_reader(file_reader)

    @staticmethod
    def _parse_file_reader(file_reader):
        hc_loci = set()
        file_reader.open()
        for line in file_reader.read_lines():
            if line.startswith("chrom\tposition\tref\tvar"):
                column_header = line
            else:
                split_line = line.split("\t")
                hc_loci.add((split_line[0], split_line[1]))
        file_reader.close()

        if not column_header:
            raise utils.JQException("Error. The hc file {} is in an incorrect"
                                    "format. Review inputs and try"
                                    "again.".format(file_reader.file_name))

        return hc_loci

    def add_tag_values(self, vcf_record):
        if (vcf_record.chrom, vcf_record.pos) not in self.hc_loci:
            if vcf_record.filter.lower() == "pass":
                vcf_record.filter = _LOW_CONFIDENCE_FILTER
            else:
                vcf_record.filter += ";" + _LOW_CONFIDENCE_FILTER
        return vcf_record

class Varscan(object):
    _HC_FILE_SUFFIX = ".Somatic.hc.fpfilter.pass"

    def __init__(self):
        self.name = "VarScan"
        self.abbr = "VS"
        self.tags = [common_tags.ReportedTag(JQ_VARSCAN_TAG),
                     common_tags.PassedTag(JQ_VARSCAN_TAG),
                     _AlleleFreqTag(),
                     _DepthTag(),
                     _SomaticTag()]
        self.meta_header = "##jacquard.normalize_varscan.sources={0},{1}\n"

    @staticmethod
    def validate_input_file(meta_headers, column_header):
        if "##source=VarScan2" not in meta_headers:
            return 0

        if VARSCAN_SOMATIC_HEADER == column_header:
            return 1
        else:
            raise utils.JQException("Unexpected VarScan VCF structure - "
                                    "missing NORMAL and TUMOR headers.")

    #TODO: (cgates) Clarify intent of this method?
    @staticmethod
    def _validate_vcf_fileset(vcf_readers):
        if len(vcf_readers) != 2:
            raise utils.JQException("VarScan directories should have exactly "
                                    "two input VCF files per patient, but "
                                    "found [{}].".format(len(vcf_readers)))

        tmp = [vcf_readers[0].file_name, vcf_readers[1].file_name]
        for i, name in enumerate(tmp):
            if "snp" in name:
                tmp[i] = "snp"
        for i, name in enumerate(tmp):
            if "indel" in name:
                tmp[i] = "indel"
        if not (tmp[0] == "snp" and tmp[1] == "indel") and not (tmp[1] == "snp"
                                                                and tmp[0] ==
                                                                "indel"):
            raise utils.JQException("Each patient in a VarScan directory "
                                    "should have a snp file and an indel file.")

        if not vcf_readers[0].column_header == vcf_readers[1].column_header:
            msg = "The column headers for VCF files [{},{}] do not match."
            raise utils.JQException(msg.format(vcf_readers[0].file_name,
                                               vcf_readers[1].file_name))

    @staticmethod
    def _validate_hc_fileset(hc_candidates):
        if len(hc_candidates) != 2:
            msg = ("VarScan directories should have exactly 2 "
                   "input somatic fpfilter files per patient, but "
                   "found [{}].")
            raise utils.JQException(msg.format(len(hc_candidates)))

        tmp = [hc_candidates[0].file_name, hc_candidates[1].file_name]
        for i, name in enumerate(tmp):
            if "snp" in name:
                tmp[i] = "snp"
        for i, name in enumerate(tmp):
            if "indel" in name:
                tmp[i] = "indel"
        if not (tmp[0] == "snp" and tmp[1] == "indel") and not (tmp[1] == "snp"
                                                                and tmp[0] ==
                                                                "indel"):
            msg = ("Each patient in a VarScan directory should "
                   "have a somatic fpfilter snp file and indel file.")
            raise utils.JQException(msg)

    def _validate_raw_input_files(self, file_readers):
        vcf_readers = []
        hc_candidates = []

        for file_reader in file_readers:
            if file_reader.file_name.lower().endswith("fpfilter.pass"):
                hc_candidates.append(file_reader)
            elif file_reader.file_name.lower().endswith(".vcf"):
                vcf_readers.append(vcf.VcfReader(file_reader))

        self._validate_vcf_fileset(vcf_readers)
        self._validate_hc_fileset(hc_candidates)

        return vcf_readers, hc_candidates

    @staticmethod
    def _parse_vcf_readers(vcf_readers, hc_keys):
        all_records = []
        metaheader_list = []
        column_header = vcf_readers[0].column_header

        for vcf_reader in vcf_readers:
            metaheader_list.extend(vcf_reader.metaheaders)
            vcf_reader.open()

#TODO: utilize the info_column dictionary in vcf
            for record in vcf_reader.vcf_records():
                if record in hc_keys:
                    record.info = record.info + ";" + JQ_VARSCAN_TAG + "HC"
                all_records.append(record)

            vcf_reader.close()

        parsed_records = [rec.asText() for rec in sorted(all_records)]

        return metaheader_list, column_header, parsed_records

    @staticmethod
    def _process_hc_files(hc_candidates):
        metaheader = None
        hc_keys = []
        for hc_file_reader in hc_candidates:
            hc_file_reader.open()

            for line in hc_file_reader.read_lines():
                split_line = line.split()

                if split_line[0] != "chrom" and split_line[0].startswith("chr"):
                    hc_key = vcf.VcfRecord(split_line[0],
                                       split_line[1],
                                       split_line[2],
                                       split_line[3])
                    hc_keys.append(hc_key)
            hc_file_reader.close()

        if len(hc_keys) > 0:
            metaheader = '##INFO=<ID=' + JQ_VARSCAN_TAG + "HC"\
                        ',Number=1,Type=Flag,Description="Jacquard '\
                        'high-confidence somatic flag for VarScan. Based on '\
                        'intersection with filtered VarScan variants.">'
        return metaheader, hc_keys

    @staticmethod
    def decorate_files(filenames, decorator):
        output_file = None
        file_name_search = "snp|indel"
        for filename in filenames:
            if not filename.lower().endswith("fpfilter.pass"):
                if re.search("("+file_name_search+")", filename):
                    prefix, suffix = re.split(file_name_search, filename)
                    output_file = os.path.basename(prefix + decorator + suffix)
                    return output_file

        raise utils.JQException("Each patient in a VarScan directory should "
                                "have a snp file and an indel file.")

    @staticmethod
    def validate_vcfs_in_directory(in_files):
        for in_file in in_files:
            if not in_file.lower().endswith("vcf"):
                if not in_file.lower().endswith("fpfilter.pass"):
                    raise utils.JQException("ERROR: Non-VCF or fpfilter file "
                                            "in directory. Check parameters "
                                            "and try again")

#TODO: Add to normalize.py
    def normalize(self, file_writer, file_readers):
        (vcf_readers,
         hc_candidates) = self._validate_raw_input_files(file_readers)

        hc_metaheader, hc_keys = self._process_hc_files(hc_candidates)
        (metaheader_list,
         column_header,
         parsed_records) = self._parse_vcf_readers(vcf_readers, hc_keys)

        if hc_metaheader is not None:
            metaheader_list.append(hc_metaheader)
        sorted_metaheader_set = sorted(set(metaheader_list))

        file_writer.open()
        for metaheader in sorted_metaheader_set:
            file_writer.write(metaheader+"\n")

        file_writer.write(column_header+"\n")

        for record in parsed_records:
            file_writer.write(record)

        file_writer.close()

    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.add_tag_values(vcf_record)
        return vcf_record.asText()

    def get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

    @staticmethod
    def _is_varscan_vcf(file_reader):
        if file_reader.file_name.endswith(".vcf"):
            vcf_reader = vcf.VcfReader(file_reader)
            return "##source=VarScan2" in vcf_reader.metaheaders
        return False

    def _is_varscan_hc_file(self, file_reader):
        return file_reader.file_name.endswith(self._HC_FILE_SUFFIX)

    def claim(self, file_readers):
        unclaimed_readers = []
        vcf_readers = []
        hc_file_readers = []
        for file_reader in file_readers:
            if self._is_varscan_vcf(file_reader):
                vcf_reader = vcf.VcfReader(file_reader)
                vcf_readers.append(vcf.RecognizedVcfReader(vcf_reader,
                                                           self))
            elif self._is_varscan_hc_file(file_reader):
                hc_file_readers.append(file_reader)
            else:
                unclaimed_readers.append(file_reader)
        return (unclaimed_readers, vcf_readers)

class _VarscanVcfReader(object):
    def __init__(self, vcf_reader, som_hc_file_reader=None):
        self._vcf_reader = vcf_reader
        self._caller = Varscan()
        self.tags = [common_tags.ReportedTag(JQ_VARSCAN_TAG),
                     common_tags.PassedTag(JQ_VARSCAN_TAG),
                     _AlleleFreqTag(),
                     _DepthTag(),
                     _SomaticTag()]

        self.som_hc_coords = None
        if som_hc_file_reader:
            self.tags.append(_HCTag(som_hc_file_reader))

    def open(self):
        return self._vcf_reader.open()

    def close(self):
        return self._vcf_reader.close()
    @property
    def metaheaders(self):
        new_metaheaders = list(self._vcf_reader.metaheaders)
        new_metaheaders.extend(self._caller.get_new_metaheaders())
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

#pylint: disable=too-few-public-methods, unused-argument
from __future__ import print_function, absolute_import
import jacquard.variant_callers.common_tags as common_tags
import jacquard.utils as utils
import jacquard.vcf as vcf
from jacquard import __version__
import re
import os

JQ_MUTECT_TAG = "JQ_MT_"

class _AlleleFreqTag(object):
    def __init__(self):
        self.metaheader = ('##FORMAT=<ID={0}AF,'
                           'Number=A,'
                           'Type=Float,'
                           #pylint: disable=line-too-long
                           'Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA)",'
                           'Source="Jacquard",'
                           'Version={1}>').format(JQ_MUTECT_TAG,
                                                  __version__)

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
    def __init__(self):
        self.metaheader = ('##FORMAT=<ID={0}DP,'
                           'Number=1,'
                           'Type=Float,'
                           #pylint: disable=line-too-long
                           'Description="Jacquard depth for MuTect (based on DP)",'
                           'Source="Jacquard",'
                           'Version={1}>').format(JQ_MUTECT_TAG,
                                                  __version__)

    @staticmethod
    def add_tag_values(vcf_record):
        if "DP" in vcf_record.format_tags:
            sample_values = {}
            for samp in vcf_record.sample_tag_values:
                sample_values[samp] = vcf_record.sample_tag_values[samp]["DP"]
            vcf_record.add_sample_tag_value(JQ_MUTECT_TAG + "DP", sample_values)

class _SomaticTag(object):
    def __init__(self):
        self.metaheader = ('##FORMAT=<ID={0}HC_SOM,'
                           'Number=1,'
                           'Type=Integer,'
                           #pylint: disable=line-too-long
                           'Description="Jacquard somatic status for MuTect: 0=non-somatic,1=somatic (based on SS FORMAT tag)",'
                           'Source="Jacquard",'
                           'Version={1}>').format(JQ_MUTECT_TAG,
                                                  __version__)

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
    _NORMAL_SAMPLE_KEY = "normal_sample_name"
    _TUMOR_SAMPLE_KEY = "tumor_sample_name"

    def __init__(self):
        self.name = "MuTect"
        self.abbr = "MT"
        self.tags = [common_tags.ReportedTag(JQ_MUTECT_TAG),
                     common_tags.PassedTag(JQ_MUTECT_TAG),
                     _AlleleFreqTag(), _DepthTag(), _SomaticTag()]
        self.file_name_search = ""

    @staticmethod
    def _get_mutect_cmd_parameters(line, mutect_dict):
        split_line = line.split(" ")

        for item in split_line:
            split_item = item.split("=")
            try:
                mutect_dict[split_item[0]] = split_item[1]
            except IndexError:
                pass

        return mutect_dict

    @staticmethod
    def decorate_files(filenames, decorator):
        output_file = None
        for file_name in filenames:
            name = re.sub(r"\.vcf$", "." + decorator + ".vcf", file_name)
            output_file = os.path.basename(name)
        return output_file

    @staticmethod
    def validate_vcfs_in_directory(in_files):
        for in_file in in_files:
            if not in_file.lower().endswith("vcf"):
                raise utils.JQException("ERROR: Non-VCF file in directory. "
                                        "Check parameters and try again")

    def normalize(self, file_writer, file_readers):
        if len(file_readers) != 1:
            raise utils.JQException(("MuTect directories should have exactly "
                                     "one input file per patient, but "
                                     "found [{}].").format(len(file_readers)))
        file_writer.open()
        for file_reader in file_readers:
            file_reader.open()

            mutect_dict = {}
            for line in file_reader.read_lines():
                if "##MuTect=" in line:
                    mutect_dict = self._get_mutect_cmd_parameters(line,
                                                                  mutect_dict)
                if "#CHROM" in line:
                    if "normal_sample_name" in mutect_dict:
                        if "tumor_sample_name" in mutect_dict:
                            line = re.sub(mutect_dict["normal_sample_name"],
                                          "NORMAL",
                                          line)
                            line = re.sub(mutect_dict["tumor_sample_name"],
                                          "TUMOR",
                                          line)
                    else:
                        raise utils.JQException("Unable to determine normal "
                                                "and tumor sample ordering "
                                                "based on MuTect metaheader.")

                file_writer.write(line)

            file_reader.close()
        file_writer.close()

    def get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

    #TODO: (cgates): Why using ints instead of boolean for this method?
    @staticmethod
    def validate_input_file(meta_headers, column_header):
        valid = 0
        for line in meta_headers:
            if "##MuTect" in line:
                valid = 1
                break
        return valid

    #TODO: (cgates): remove this when translate is complete
    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.add_tag_values(vcf_record)
        return vcf_record.asText()

    @staticmethod
    def _is_mutect_vcf(file_reader):
        if file_reader.file_name.endswith(".vcf"):
            vcf_reader = vcf.VcfReader(file_reader)
            for metaheader in vcf_reader.metaheaders:
                if metaheader.startswith("##MuTect"):
                    return True
        return False

    def claim(self, file_readers):
        unclaimed_readers = []
        vcf_readers = []
        for file_reader in file_readers:
            if self._is_mutect_vcf(file_reader):
                vcf_reader = vcf.VcfReader(file_reader)
                vcf_readers.append(_MutectVcfReader(vcf_reader))
            else:
                unclaimed_readers.append(file_reader)
        return (unclaimed_readers, vcf_readers)

    def _add_tags(self, vcf_record):
        for tag in self.tags:
            tag.add_tag_values(vcf_record)
        return vcf_record

    @staticmethod
    def _build_mutect_dict(metaheaders):
        mutect_dict = {}
        for metaheader in metaheaders:
            if metaheader.startswith("##MuTect="):
                split_line = metaheader.strip('"').split(" ")
                for item in split_line:
                    split_item = item.split("=")
                    try:
                        mutect_dict[split_item[0]] = split_item[1]
                    except IndexError:
                        pass

        return mutect_dict

    def _get_new_column_header(self, vcf_reader):
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

    def vcf_records(self):
        for vcf_record in self._vcf_reader.vcf_records():
            yield self._add_tags(vcf_record)

    def _add_tags(self, vcf_record):
        for tag in self.tags:
            tag.add_tag_values(vcf_record)
        return vcf_record

    def add_tag_class(self, tag_classes):
        self.tags.extend(tag_classes)


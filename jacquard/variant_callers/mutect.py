#pylint: disable=too-few-public-methods, unused-argument
from __future__ import print_function, absolute_import
import jacquard.variant_callers.common_tags as common_tags
import jacquard.utils as utils
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

    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.add_tag_values(vcf_record)
        return vcf_record.asText()

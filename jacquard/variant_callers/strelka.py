from __future__ import print_function, absolute_import
from jacquard.vcf import VcfReader
import jacquard.utils as utils
import os
import re

JQ_STRELKA_TAG = "JQ_SK_"

#pylint: disable=too-few-public-methods
class _AlleleFreqTag(object):
    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}AF,'
                           'Number=A,'
                           'Type=Float,'
                           'Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth. Uses (TIR tier 2)/DP2 if available, otherwise uses (ACGT tier2 depth) / DP2)",'
                           'Source="Jacquard",'
                           'Version={1}>').format(JQ_STRELKA_TAG,
                                                  utils.__version__)

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

            rounded_af = self._round_two_digits(str(freq))
            capped_af = min(rounded_af, "1.00")
            afs.append(capped_af)

        return afs

    def _get_indelallelefreq_per_sample(self, vcf_record, sample):
        afs = []
        indel_depths = vcf_record.sample_tag_values[sample]["TIR"]
        numerator = float(indel_depths.split(",")[1])
        denominator = float(vcf_record.sample_tag_values[sample]["DP2"])
        freq = numerator/denominator if denominator != 0 else 0.0

        rounded_af = self._round_two_digits(str(freq))
        capped_af = min(rounded_af, "1.00")
        afs.append(capped_af)

        return afs

    def format(self, vcf_record):
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
    def _round_two_digits(value):
        if len(value.split(".")[1]) <= 2:
            return value
        else:
            return str(round(100 * float(value))/100)

class _DepthTag(object):
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
                           'Description="Jacquard depth for Strelka (uses DP2 if available, otherwise uses ACGT tier2 depth)",'
                           'Source="Jacquard",'
                           'Version={1}>').format(JQ_STRELKA_TAG,
                                                  utils.__version__)

    def format(self, vcf_record):
        if vcf_record.format_tags.isdisjoint(_DepthTag.REQUIRED_TAGS):
            return
        sample_values = {}
        for sample in vcf_record.sample_tag_values:
            sample_tags = vcf_record.sample_tag_values[sample]
            sample_values[sample] = _DepthTag._get_tier2_base_depth(sample_tags)
        vcf_record.add_sample_tag_value(JQ_STRELKA_TAG + "DP", sample_values)

class _SomaticTag(object):
    def __init__(self):
        #pylint: disable=line-too-long
        self.metaheader = ('##FORMAT=<ID={0}HC_SOM,'
                           'Number=1,Type=Integer,'
                           'Description="Jacquard somatic status for Strelka: 0=non-somatic,1=somatic (based on PASS in FILTER column)",'
                           'Source="Jacquard",'
                           'Version={1}>').format(JQ_STRELKA_TAG,
                                                  utils.__version__)

    def format(self, vcf_record):
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
    def __init__(self):
        self.name = "Strelka"
        self.tags = [_AlleleFreqTag(), _DepthTag(), _SomaticTag()]
        self.meta_header = "##jacquard.normalize_strelka.sources={0},{1}\n"

    @staticmethod
    def _validate_raw_input_files(file_readers):
        if len(file_readers) != 2:
            raise utils.JQException("Strelka directories should have exactly "
                                    "two input files per patient, but found "
                                    "[{}].".format(len(file_readers)))

        tmp = [file_readers[0].file_name, file_readers[1].file_name]
        for i, name in enumerate(tmp):
            if "snvs" in name:
                tmp[i] = "snvs"
        for i, name in enumerate(tmp):
            if "indels" in name:
                tmp[i] = "indels"
        if not (tmp[0] == "snvs" and tmp[1] == "indels") and not (tmp[1] ==
                                                                  "snvs" and
                                                                  tmp[0] ==
                                                                  "indels"):
            raise utils.JQException("Each patient in a Strelka directory "
                                    "should have a snvs file and an indels "
                                    "file.")

        vcf_readers = [VcfReader(file_readers[0]), VcfReader(file_readers[1])]
        if not vcf_readers[0].column_header == vcf_readers[1].column_header:
            raise utils.JQException("The column headers for VCF files [{},{}] "
                                    "do not match."\
                .format(vcf_readers[0].file_name, vcf_readers[1].file_name))
        return vcf_readers

    @staticmethod
    def _parse_vcf_readers(vcf_readers):
        all_records = []
        metaheader_list = []
        column_header = vcf_readers[0].column_header

        for vcf_reader in vcf_readers:
            metaheader_list.extend(vcf_reader.metaheaders)
            vcf_reader.open()

            for record in vcf_reader.vcf_records():
                all_records.append(record)
            vcf_reader.close()

        parsed_records = [rec.asText() for rec in sorted(all_records)]

        return metaheader_list, column_header, parsed_records

    def decorate_files(self, filenames, decorator):
        output_file = None
        file_name_search = "snvs|indels"
        for filename in filenames:
            if re.search("("+file_name_search+")", filename):
                prefix, suffix = re.split(file_name_search, filename)
                output_file = os.path.basename(prefix+decorator+suffix)
                return output_file

        raise utils.JQException("Each patient in a Strelka directory should "
                                "have a snvs file and an indels file.")

    def validate_vcfs_in_directory(self, in_files):
        for in_file in in_files:
            if not in_file.lower().endswith("vcf"):
                raise utils.JQException("ERROR: Non-VCF file in directory. "
                                        "Check parameters and try again")

    def normalize(self, file_writer, file_readers):
        vcf_readers = self._validate_raw_input_files(file_readers)
        metaheader_list, column_header, parsed_records = self._parse_vcf_readers(vcf_readers)
        sorted_metaheader_set = sorted(set(metaheader_list))

        file_writer.open()

        for metaheader in sorted_metaheader_set:
            file_writer.write(metaheader+"\n")
        file_writer.write(column_header+"\n")

        for record in parsed_records:
            file_writer.write(record)

        file_writer.close()

    def get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

    def validate_input_file(self, meta_headers, column_header):
        return "##source=strelka" in meta_headers

    def validate_record(self,vcfRecord):
        return True

    def final_steps(self, hc_candidates, merge_candidates, output_dir):
        output_file_count = len(merge_candidates.keys())
        print ("Wrote [{0}] VCF files to [{1}]").format(output_file_count,
                                                        output_dir)
        return merge_candidates

    def handle_hc_files(self, in_file, out_dir, hc_candidates):
        return hc_candidates

    def validate_file_set(self, all_keys):
        pass

    def add_tags(self, vcf_record):
        for tag in self.tags:
            tag.format(vcf_record)
        return vcf_record.asText()

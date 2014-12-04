import jacquard.utils as utils
from jacquard.vcf import VcfReader
import re
from jacquard.utils import JQException
import os

JQ_STRELKA_TAG = "JQ_SK_"

class _AlleleFreqTag(object):
    def __init__(self):
        self.metaheader = ('##FORMAT=<ID={0}AF,'
                           'Number=A,'
                           'Type=Float,'
                           'Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth. Uses (TIR tier 2)/DP2 if available, otherwise uses (ACGT tier2 depth) / DP2)",'
                           'Source="Jacquard",'
                           'Version={1}>').format(JQ_STRELKA_TAG, 
                                                  utils.__version__)

    def _get_tier2_base_depth(self, sample_format_dict, alt_allele):
        numerator = float(sample_format_dict[alt_allele + "U"].split(",")[1])
        tags = ["AU", "CU", "TU", "GU"]
        depth = 0

        for tag in tags:
            depth += float(sample_format_dict[tag].split(",")[1])
        af = numerator/depth if depth != 0 else 0.0

        return af

    def _get_SNPallelefreq_per_sample(self, vcfRecord, key):
        afs = []
        split_alt = vcfRecord.alt.split(",")

        for alt_allele in split_alt:
            sample_format_dict = vcfRecord.sample_dict[key]
            af = self._get_tier2_base_depth(sample_format_dict, alt_allele)

            rounded_af = self._round_two_digits(str(af))
            capped_af = min(rounded_af, "1.00")
            afs.append(capped_af)

        return afs

    def _get_indelallelefreq_per_sample(self, vcfRecord, key):
        afs = []
        numerator = float(vcfRecord.sample_dict[key]["TIR"].split(",")[1])
        denominator = float(vcfRecord.sample_dict[key]["DP2"])
        af = numerator/denominator if denominator != 0 else 0.0

        rounded_af = self._round_two_digits(str(af))
        capped_af = min(rounded_af, "1.00")
        afs.append(capped_af)

        return afs 

    def format(self, vcfRecord):
        sample_values = {}
        if vcfRecord.alt == ".":
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = "."
        else:
            for key in vcfRecord.sample_dict.keys():
                if "AU" in vcfRecord.format_set:#if it's an snp
                    afs = self._get_SNPallelefreq_per_sample(vcfRecord,key)
                    sample_values[key] = ",".join(afs)
                elif "TIR" in vcfRecord.format_set: #if it's an indel
                    afs = self._get_indelallelefreq_per_sample(vcfRecord,key)
                    sample_values[key] = ",".join(afs)
                else:
                    continue

        if sample_values:
            vcfRecord.insert_format_field(JQ_STRELKA_TAG + "AF" ,sample_values)

    def _round_two_digits(self, value):
        if len(value.split(".")[1]) <= 2:
            return value
        else:
            return str(round(100 * float(value))/100)

class _DepthTag(object):
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}DP,Number=1,Type=Float,Description="Jacquard depth for Strelka (uses DP2 if available, otherwise uses ACGT tier2 depth)",Source="Jacquard",Version={1}>'.format(JQ_STRELKA_TAG, utils.__version__)

    def _get_tier2_base_depth(self, sample_format_dict):
        tags = ["AU", "CU", "TU", "GU"]
        depth = 0
        for tag in tags:
            depth += int(sample_format_dict[tag].split(",")[1])
        return depth

    def format(self, vcfRecord):
        if "DP2" not in vcfRecord.format_set and "AU" not in vcfRecord.format_set:
            return

        sample_values = {}
        if "DP2" in vcfRecord.format_set:
            for key in vcfRecord.sample_dict:
                sample_values[key] = vcfRecord.sample_dict[key]["DP2"]
        elif "AU" in vcfRecord.format_set:
            for key in vcfRecord.sample_dict:
                sample_format_dict = vcfRecord.sample_dict[key]
                sample_values[key] = self._get_tier2_base_depth(sample_format_dict)

        vcfRecord.insert_format_field(JQ_STRELKA_TAG + "DP", sample_values)

class _SomaticTag(object):
    #TODO: cgates :Pull tag metaheaders to resource bundle?
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}HC_SOM,Number=1,Type=Integer,Description="Jacquard somatic status for Strelka: 0=non-somatic,1=somatic (based on PASS in FILTER column)",Source="Jacquard",Version={1}>'.format(JQ_STRELKA_TAG, utils.__version__)

    # pylint: disable=W0613,R0201
    def format(self, vcfRecord):
        strelka_tag = JQ_STRELKA_TAG + "HC_SOM"
        sample_values = {}
        if vcfRecord.filter == "PASS":
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = self._somatic_status(key)
        else:
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = "0"

        vcfRecord.insert_format_field(strelka_tag,sample_values)

    def _somatic_status(self,count):
        if count == 0: #it's NORMAL
            return "0"
        else: #it's TUMOR
            return "1"

class Strelka(object):
    def __init__(self):
        self.name = "Strelka"
        self.good = True
        self.tags = [_AlleleFreqTag(),_DepthTag(),_SomaticTag()]
        self.meta_header = "##jacquard.normalize_strelka.sources={0},{1}\n"

#TODO: Add to normalize.py.
    def _validate_raw_input_files(self, file_readers):
        if len(file_readers) != 2:
            raise utils.JQException("Strelka directories should have exactly "
                                    "two input files per patient, but found "
                                    "[{}].".format(len(file_readers)))

        tmp = [file_readers[0].file_name,file_readers[1].file_name]
        for i,name in enumerate(tmp):
            if "snvs" in name:
                tmp[i] = "snvs"
        for i,name in enumerate(tmp):
            if "indels" in name:
                tmp[i] = "indels"
        if not (tmp[0] == "snvs" and tmp[1] == "indels") and not (tmp[1] ==
                                                                  "snvs" and
                                                                  tmp[0] ==
                                                                  "indels"):
            raise utils.JQException("Each patient in a Strelka directory should "
                                    "have a snvs file and an indels file.")

        vcf_readers = [VcfReader(file_readers[0]),VcfReader(file_readers[1])]
        if not vcf_readers[0].column_header == vcf_readers[1].column_header:
            raise utils.JQException("The column headers for VCF files [{},{}] "
                                    "do not match."\
                .format(vcf_readers[0].file_name,vcf_readers[1].file_name))
        return vcf_readers

    def _parse_vcf_readers(self,vcf_readers):
        all_records = []
        metaheader_list = []
        column_header = vcf_readers[0].column_header

        for vcf_reader in vcf_readers:
            metaheader_list.extend(vcf_reader.metaheaders)
            vcf_reader.open()

            for record in vcf_reader.vcf_records():
                all_records.append(record.asText())
            vcf_reader.close()

        parsed_records = utils.sort_data(all_records)

        return metaheader_list, column_header, parsed_records

    def decorate_files(self, filenames, decorator):
        output_file = None
        file_name_search = "snvs|indels"
        for filename in filenames:
            if re.search("("+file_name_search+")", filename):
                prefix,suffix = re.split(file_name_search,filename)
                output_file = os.path.basename(prefix+decorator+suffix)
                return output_file

        raise utils.JQException("Each patient in a Strelka directory should "
                                "have a snvs file and an indels file.")

    def validate_vcfs_in_directory(self, in_files):
        for in_file in in_files:
            if not in_file.lower().endswith("vcf"):
                raise utils.JQException("ERROR: Non-VCF file in directory. "
                                        "Check parameters and try again")

#TODO: Add to normalize.py.        
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
        print "Wrote [{0}] VCF files to [{1}]". \
            format(len(merge_candidates.keys()), output_dir)
        return merge_candidates

    def handle_hc_files(self, in_file, out_dir, hc_candidates):
        return hc_candidates

    def validate_file_set(self, all_keys):
        pass

    def add_tags(self,vcfRecord):
        for tag in self.tags:
            tag.format(vcfRecord)
        return vcfRecord.asText()

from collections import defaultdict
import os
import re
import jacquard.utils as utils
from jacquard.vcf import VcfReader


_VARSCAN_SOMATIC_HEADER = "#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR".replace("|","\t")
JQ_VARSCAN_TAG = "JQ_VS_"

class AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}AF,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ)",Source="Jacquard",Version={1}>'.format(JQ_VARSCAN_TAG, utils.__version__)

    def format(self, vcfRecord):
        sample_values = {}
        if "FREQ" in vcfRecord.format_set:
            for key in vcfRecord.sample_dict.keys():
                freq = vcfRecord.sample_dict[key]["FREQ"].split(",")
                sample_values[key] = self.roundTwoDigits(freq)
            vcfRecord.insert_format_field(JQ_VARSCAN_TAG + "AF",sample_values)

    def roundTwoDigits(self, value): 
        new_values = []
        for val in value:
            new_val = str(float(val.strip("%"))/100)
            if len(new_val.split(".")[1]) <= 2:
                new_values.append(new_val)
            else:
                new_values.append(str(round(100 * float(new_val))/100))
        return ",".join(new_values) 

class DepthTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}DP,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP)",Source="Jacquard",Version={1}>'.format(JQ_VARSCAN_TAG, utils.__version__)

    def format(self, vcfRecord):
        if "DP" in vcfRecord.format_set:
            sample_values = {}
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = vcfRecord.sample_dict[key]["DP"]
            vcfRecord.insert_format_field(JQ_VARSCAN_TAG + "DP",sample_values)

class SomaticTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}HC_SOM,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1=somatic (based on SOMATIC info tag and if sample is TUMOR)",Source="Jacquard",Version={1}>'.format(JQ_VARSCAN_TAG, utils.__version__)

    def format(self, vcfRecord):
        info_array = vcfRecord.info.split(";")
        varscan_tag = JQ_VARSCAN_TAG + "HC_SOM"
        sample_values = {}
        if "SS=2" in info_array and JQ_VARSCAN_TAG + "HC" in info_array:
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = self.somatic_status(key)
        else:
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = "0"

        vcfRecord.insert_format_field(varscan_tag,sample_values)

    def somatic_status(self, count):
        if count == 0: #it's NORMAL
            return "0"
        else: #it's TUMOR
            return "1"

class Varscan():
    def __init__(self):
        self.name = "VarScan"
        self.good = True
        self.tags = [AlleleFreqTag(),DepthTag(),SomaticTag()]
        self.meta_header = "##jacquard.normalize_varscan.sources={0},{1}\n"

    def validate_input_file(self, meta_headers, column_header):
        valid = 0
        for line in meta_headers:
            if "##source=VarScan2" in line:
                valid = 1
        if not valid:
            return 0
        
        if _VARSCAN_SOMATIC_HEADER == column_header:
            return 1
        else:
            raise utils.JQException("Unexpected VarScan VCF structure - "
                                    "missing NORMAL and TUMOR headers.")

    def _validate_vcf_fileset(self, vcf_readers):
        if len(vcf_readers) != 2:
            raise utils.JQException("VarScan directories should have exactly "
                                    "two input VCF files per patient, but "
                                    "found [{}].".format(len(vcf_readers)))

        tmp = [vcf_readers[0].file_name, vcf_readers[1].file_name]
        for i,name in enumerate(tmp):
            if "snp" in name:
                tmp[i] = "snp"
        for i,name in enumerate(tmp):
            if "indel" in name:
                tmp[i] = "indel"
        if not (tmp[0] == "snp" and tmp[1] == "indel") and not (tmp[1] == "snp"
                                                                and tmp[0] ==
                                                                "indel"):
            raise utils.JQException("Each patient in a VarScan directory "
                                    "should have a snp file and an indel file.")

        if not vcf_readers[0].column_header == vcf_readers[1].column_header:
            raise utils.JQException("The column headers for VCF files [{},{}] "
                                    "do not match."\
                .format(vcf_readers[0].file_name,vcf_readers[1].file_name))

    def _validate_hc_fileset(self, hc_candidates):
        if len(hc_candidates) != 2:
            raise utils.JQException("VarScan directories should have exactly 2 "
                                    "input somatic fpfilter files per patient, but "
                                    "found [{}].".format(len(hc_candidates)))

        tmp = [hc_candidates[0].file_name,hc_candidates[1].file_name]
        for i,name in enumerate(tmp):
            if "snp" in name:
                tmp[i] = "snp"
        for i,name in enumerate(tmp):
            if "indel" in name:
                tmp[i] = "indel"
        if not (tmp[0] == "snp" and tmp[1] == "indel") and not (tmp[1] == "snp"
                                                                and tmp[0] ==
                                                                "indel"):
            raise utils.JQException("Each patient in a VarScan directory should "
                                    "have a somatic fpfilter snp file and indel file.")

        pass

    def _validate_raw_input_files(self, file_readers):
        vcf_readers = []
        hc_candidates = []

        for file_reader in file_readers:
            if file_reader.file_name.lower().endswith("fpfilter.pass"):
                hc_candidates.append(file_reader)
            elif file_reader.file_name.lower().endswith(".vcf"):
                vcf_readers.append(VcfReader(file_reader))

        self._validate_vcf_fileset(vcf_readers)
        self._validate_hc_fileset(hc_candidates)

        return vcf_readers, hc_candidates

    def _parse_vcf_readers(self,vcf_readers,hc_keys):
        all_records = []
        metaheader_list = []
        column_header = vcf_readers[0].column_header

        for vcf_reader in vcf_readers:
            metaheader_list.extend(vcf_reader.metaheaders)
            vcf_reader.open()

            for record in vcf_reader.vcf_records():
                if record.key in hc_keys:
                    record.info = record.info + ";" + JQ_VARSCAN_TAG + "HC"
                all_records.append(record.asText())

            vcf_reader.close()

        parsed_records = utils.sort_data(all_records)
        return metaheader_list, column_header, parsed_records
    
    def _process_hc_files(self, hc_candidates):
        metaheader = None
        hc_keys = []
        for hc_file_reader in hc_candidates:
            hc_file_reader.open()

            for line in hc_file_reader.read_lines():
                split_line = line.split()

                if split_line[0] != "chrom" and split_line[0].startswith("chr"):
                    hc_key = split_line[0] + "_" + split_line[1] + "_" + split_line[2] + "_" + split_line[3]
                    hc_keys.append(hc_key)
            hc_file_reader.close()

        if len(hc_keys)>0:
            metaheader = '##INFO=<ID=' + JQ_VARSCAN_TAG + "HC"\
                        ',Number=1,Type=Flag,Description="Jaquard '\
                        'high-confidence somatic flag for VarScan. Based on '\
                        'intersection with filtered VarScan variants.">'
        return metaheader, hc_keys

    def decorate_files(self, filenames, decorator):
        output_file = None
        file_name_search = "snp|indel"
        for filename in filenames:
            if not filename.lower().endswith("fpfilter.pass"):
                if re.search("("+file_name_search+")", filename):
                    prefix,suffix = re.split(file_name_search,filename)
                    output_file = os.path.basename(prefix+decorator+suffix)
                    return output_file

        raise utils.JQException("Each patient in a VarScan directory should "
                                "have a snp file and an indel file.")
    
    def validate_vcfs_in_directory(self, in_files):
        for in_file in in_files:
            if not in_file.lower().endswith("vcf") and not in_file.lower().endswith("fpfilter.pass"):
                raise utils.JQException("ERROR: Non-VCF or fpfilter file "
                                        "in directory. Check parameters and "
                                        "try again")

#TODO: Add to normalize.py.        
    def normalize(self, file_writer, file_readers):
        vcf_readers, hc_candidates = self._validate_raw_input_files(file_readers)

        hc_metaheader, hc_keys = self._process_hc_files(hc_candidates)
        metaheader_list, column_header, parsed_records = self._parse_vcf_readers(vcf_readers, hc_keys)

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

    def add_tags(self,vcfRecord):
        for tag in self.tags:
            tag.format(vcfRecord)
        return vcfRecord.asText()

    def get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

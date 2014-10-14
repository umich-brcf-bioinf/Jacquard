from collections import defaultdict
import os
import re
import jacquard.utils as utils
from jacquard.vcf import VcfReader


_VARSCAN_SOMATIC_HEADER = "#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|NORMAL|TUMOR".replace("|","\t")
_JQ_VARSCAN_HC_INFO_FIELD = "JQ_HC_VS"

class AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}VS,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ)",Source="Jacquard",Version={1}>'.format(utils.jq_af_tag, utils.__version__)

    def format(self, vcfRecord):
        sample_values = {}
        if "FREQ" in vcfRecord.format_set:
            for key in vcfRecord.sample_dict.keys():
                freq = vcfRecord.sample_dict[key]["FREQ"].split(",")
                sample_values[key] = self.roundTwoDigits(freq)
            vcfRecord.insert_format_field("JQ_AF_VS",sample_values)
                
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

        self.metaheader = '##FORMAT=<ID={0}VS,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP)",Source="Jacquard",Version={1}>'.format(utils.jq_dp_tag, utils.__version__)

    def format(self, vcfRecord):
        if "DP" in vcfRecord.format_set:
            sample_values = {}
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = vcfRecord.sample_dict[key]["DP"]
            vcfRecord.insert_format_field("JQ_DP_VS",sample_values)
    
class SomaticTag():
    def __init__(self):

        self.metaheader = '##FORMAT=<ID={0}VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1=somatic (based on SOMATIC info tag and if sample is TUMOR)",Source="Jacquard",Version={1}>'.format(utils.jq_somatic_tag, utils.__version__)

    def format(self, vcfRecord):
        info_array = vcfRecord.info.split(";")
        varscan_tag = utils.jq_somatic_tag + "VS"
        sample_values = {}
        if "SS=2" in info_array and _JQ_VARSCAN_HC_INFO_FIELD in info_array:
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
        self.file_name_search = "snp|indel"
        
    def validate_input_file(self, meta_headers, column_header):
        if "##source=VarScan2" not in meta_headers:
            return 0

        if _VARSCAN_SOMATIC_HEADER == column_header:
            return 1
        else:
            raise utils.JQException("Unexpected VarScan VCF structure - missing NORMAL and TUMOR headers.")

    def _validate_vcf_fileset(self,vcf_readers):
        if len(vcf_readers) != 2:
            raise utils.JQException("ERROR: VarScan directories should have exactly two input VCF files per patient, but found [{}].".format(len(vcf_readers)))

        tmp = [vcf_readers[0].file_name,vcf_readers[1].file_name]
        for i,name in enumerate(tmp):
            if "snp" in name:
                tmp[i] = "snp"
        for i,name in enumerate(tmp):
            if "indel" in name:
                tmp[i] = "indel"
        if not (tmp[0] == "snp" and tmp[1] == "indel") and not (tmp[1] == "snp" and tmp[0] == "indel"): 
            raise utils.JQException("ERROR: Each patient in a VarScan directory should have a snp file and an indel file.")

        if not vcf_readers[0].column_header == vcf_readers[1].column_header:
            raise utils.JQException("ERROR: The column headers for VCF files [{},{}] do not match."\
                .format(vcf_readers[0].file_name,vcf_readers[1].file_name))
    
    def _validate_hc_fileset(self,hc_candidates):
        if len(hc_candidates) != 6:
            raise utils.JQException("ERROR: VarScan directories should have exactly 6 input HC files per patient, but found [{}].".format(len(hc_candidates)))
        pass
        
    def _validate_raw_input_files(self, file_readers):
        vcf_readers = []
        hc_candidates = []
        for file_reader in file_readers:
            fname, ext = os.path.splitext(file_reader.file_name)
            if ext == ".hc":
                hc_candidates.append(file_reader)
            elif ext == ".vcf":
                vcf_readers.append(VcfReader(file_reader))
                
        self._validate_vcf_fileset(vcf_readers)
        self._validate_hc_fileset(hc_candidates)
              
        return vcf_readers,hc_candidates
        
    def _write_vcf_records(self,vcf_readers):
        all_records = []
        for vcf_reader in vcf_readers:
            for record in vcf_reader.vcf_records():
                all_records.append(record.asText())
        parsed_records = utils.sort_data(all_records)
        return parsed_records
        
#     def _identify_hc_records(self, file_readers):
#         hc_candidates = []
#         for file_reader in file_readers:
#             fname, ext = os.path.splitext(file_reader.file_name)
#             if ext == ".hc":
#                 hc_candidates.append()
                
#TODO: Add to normalize.py.        
    def normalize(self, file_writer, file_readers):
        vcf_readers,hc_candidates = self._validate_raw_input_files(file_readers)
        metaheader_list = []
        column_header = None
        for i,vcf_reader in enumerate(vcf_readers):
            if i==0:
                column_header = vcf_reader.column_header
            metaheader_list.extend(vcf_reader.metaheaders)
        sorted_metaheader_set = sorted(set(metaheader_list))
        file_writer.open()
        for metaheader in sorted_metaheader_set:
            file_writer.write(metaheader+"\n")
        file_writer.write(column_header+"\n")
        parsed_records = self._write_vcf_records(vcf_readers)
        for record in parsed_records:
            file_writer.write(record)
        file_writer.close()
        
    def identify_hc_variants(self,hc_candidates):
        hc_variants = {}
        for key, vals in hc_candidates.items():
            for file in vals:
                f = open(file, "r")
                for line in f:
                    split_line = line.split("\t")
                    if line.startswith("chrom"):
                        continue
                    else:
                        hc_key = "^".join([split_line[0], split_line[1], split_line[2], split_line[3]])
                        hc_variants[hc_key] = 1
                
        return hc_variants

    def mark_hc_variants(self,hc_variants, merge_candidates, output_dir):
        marked_as_hc = []
    
        for key, vals in merge_candidates.items():
            new_lines = []
            headers = []
            f = open(key, "r")
            for line in f:
                split_line = line.split("\t")
                if line.startswith('"'):
                    line = line.replace("\t", "")
                    line = line[1:-2] + "\n"
                if line.startswith("#"):
                    headers.append(line)
                else:
                    merge_key = "^".join([split_line[0], str(split_line[1]), split_line[3], split_line[4]])
                    if merge_key in hc_variants:
                        if _JQ_VARSCAN_HC_INFO_FIELD not in split_line[7]:
                            split_line[7] += ";"+_JQ_VARSCAN_HC_INFO_FIELD
                        marked_as_hc.append(merge_key)
                    new_line = "\t".join(split_line)
                    new_lines.append(new_line)
            f.close()
            
            sorted_headers = utils.sort_headers(headers)
            self.write_to_merged_file(new_lines, sorted_headers, key)
        
        print "Wrote [{0}] VCF files to [{1}]".format(len(merge_candidates.keys()), output_dir)
        
        return marked_as_hc
    
    def write_to_merged_file(self, new_lines, headers, key):
        sorted_variants = utils.sort_data(new_lines)
        
        writer = open(key, "w")
        utils.write_output(writer, headers, sorted_variants)
        writer.close()

    def final_steps(self, hc_candidates, merge_candidates, output_dir):
        hc_variants = self.identify_hc_variants(hc_candidates)
        marked_as_hc = self.mark_hc_variants(hc_variants, merge_candidates, output_dir)

        return marked_as_hc
    
    def handle_hc_files(self, in_file, out_dir, hc_candidates):
        merged_fname = re.sub("snp|indel", "merged", os.path.join(out_dir, os.path.basename(in_file)))
        hc_candidates[merged_fname].append(in_file)
        
        return hc_candidates
        
    def validate_file_set(self, all_keys):
        sample_files = defaultdict(list)
        for key in all_keys:
            prefix = key.split("merged")[0]
            suffix = key.split("merged")[1]
            sample_files[prefix.strip("_")].append(suffix)
    
        required_vals = [".Germline.hc", ".LOH.hc", ".Somatic.hc", ".vcf"]
        missing = 0
        added = 0
        for key, val in sample_files.items():
            missing_files, missing = self.check_for_missing(required_vals, val, key, missing)
            added_files, added = self.check_for_unknown(required_vals, val, key, added)
            
        if missing == 1:
            print "ERROR: Some required files were missing. Review input directory and try again"
            exit(1)
        if added == 1:
            print "WARNING: Some samples had unknown .hc files"
            
        return sample_files

    def check_for_missing(self, required_vals, val, key, missing):
        missing_files = []
        for item in required_vals:
            if item not in val:
                missing_files.append(item)
                missing = 1
        if missing_files != []:
            print "ERROR: [{0}] is missing required files {1}".format(key.strip("_"), missing_files)
         
        return missing_files, missing
     
    def check_for_unknown(self, required_vals, val, key, added):
        added_files = []
        for thing in val:
            if thing not in required_vals:
                added_files.append(thing)
                added = 1
        if added_files != []:
            print "WARNING: [{0}] has unknown .hc files {1}".format(key.strip("_"), added_files)
         
        return added_files, added
    
    def add_tags(self,vcfRecord):
        for tag in self.tags:
            tag.format(vcfRecord)
        return vcfRecord.asText()
    
    def get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

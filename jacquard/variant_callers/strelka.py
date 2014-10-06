import jacquard.jacquard_utils as jacquard_utils

class _AlleleFreqTag(object):
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}SK,Number=A,Type=Float,Description="Jacquard allele frequency for Strelka: Decimal allele frequency rounded to 2 digits (based on alt_depth/total_depth. Uses TAR if available, otherwise uses uses DP2 if available, otherwise uses ACGT tier2 depth)",Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_af_tag, jacquard_utils.__version__)

    def format(self, vcfRecord):
        afs = []
        if vcfRecord.alt == ".":
            afs = ["."]
        else:
            split_alt = vcfRecord.alt.split(",")
            for alt_allele in split_alt:
                if "AU" in vcfRecord.format_set: #if it's an snv
                    for key in vcfRecord.sample_dict.keys():
                        numerator = float(vcfRecord.sample_dict[key][alt_allele + "U"].split(",")[1])
                        tags = ["AU", "CU", "TU", "GU"]
                        denominator = 0
                        for tag in tags:
                            denominator += float(vcfRecord.sample_dict[key][tag].split(",")[1])
                        af = numerator/denominator if denominator != 0 else 0.0
                elif "TAR" in vcfRecord.format_set: #if it's an indel
                    for key in vcfRecord.sample_dict.keys():
                        numerator = float(vcfRecord.sample_dict[key]["TAR"].split(",")[1])
                        denominator = float(vcfRecord.sample_dict[key]["DP2"])
                        af = numerator/denominator if denominator != 0 else 0.0
                else:
                    continue

                rounded_af = self._round_two_digits(str(af))
                capped_af = min(rounded_af, "1.00")
                afs.append(capped_af)

        if afs != []:
            for key in vcfRecord.sample_dict.keys():
                vcfRecord.sample_dict[key]["JQ_AF_SK"] = ",".join(afs)

    def _round_two_digits(self, value):
        if len(value.split(".")[1]) <= 2:
            return value
        else:
            return str(round(100 * float(value))/100)

class _DepthTag(object):
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}SK,Number=1,Type=Float,Description="Jacquard depth for Strelka (uses DP2 if available, otherwise uses ACGT tier2 depth),Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_dp_tag, jacquard_utils.__version__)

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

        vcfRecord.insert_format_field("JQ_DP_SK",sample_values)
        


class _SomaticTag(object):
    #TODO: cgates :Pull tag metaheaders to resource bundle?
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}SK,Number=1,Type=Integer,Description="Jacquard somatic status for Strelka: 0=non-somatic,1= somatic (based on PASS in FILTER column),Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_somatic_tag, jacquard_utils.__version__)

    # pylint: disable=W0613,R0201
    #TODO: cgates : Refactor params to record_object?
    def format(self, vcfRecord):
        strelka_tag = jacquard_utils.jq_somatic_tag + "SK"
        if vcfRecord.filter_field == "PASS":
            for key in vcfRecord.sample_dict.keys():
                vcfRecord.sample_dict[key][strelka_tag] = self._somatic_status(key)
        else:
            for key in vcfRecord.sample_dict.keys():
                vcfRecord.sample_dict[key][strelka_tag] = "0"

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
        self.file_name_search = "snvs|indels"

    def validate_input_file(self, header):
        valid = 0
        for line in header:
            if line.startswith("##source=strelka"):
                valid = 1
            elif line.startswith("##"):
                continue
            else:
                break
        return (valid)
    
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
    
    def update_metaheader(self,metaheader):
        for tag in self.tags:
            metaheader += tag.metaheader
        return metaheader
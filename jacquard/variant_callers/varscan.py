class Varscan():
    def __init__(self):
        self.name = "VarScan"
        self.meta_header = "##jacquard.normalize_varscan.sources={0},{1}\n"
        self.file_name_search = "snp|indel"
        
    def validate_input_file(self, input_file):
        valid = 0
        for line in input_file:
            if line.startswith("##source=VarScan2"):
                valid = 1
            elif line.startswith("##"):
                continue
            else:
                break
        return (self.name, valid)
    
    def final_steps(self, hc_candidates, merge_candidates, output_dir):
        hc_variants = identify_hc_variants(hc_candidates)
        marked_as_hc = mark_hc_variants(hc_variants, merge_candidates, output_dir)

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
            missing_files, missing = check_for_missing(required_vals, val, key, missing)
            added_files, added = check_for_unknown(required_vals, val, key, added)
            
        if missing == 1:
            print "ERROR: Some required files were missing. Review input directory and try again"
            exit(1)
        if added == 1:
            print "WARNING: Some samples had unknown .hc files"
            
        return sample_files
    
    class Varscan_AlleleFreqTag():
        def __init__(self):
            self.metaheader = '##FORMAT=<ID=JQ_AF_VS,Number=A,Type=Float,Description="Jacquard allele frequency for VarScan: Decimal allele frequency rounded to 2 digits (based on FREQ),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)
    
        def format(self, alt, filter, info_string, format_dict, count):
            if "FREQ" in format_dict.keys():
                freq = format_dict["FREQ"].split(",")
                format_dict["JQ_AF_VS"] = self.roundTwoDigits(freq)
                
            return format_dict
    
        def roundTwoDigits(self, value): 
            new_values = []
            for val in value:
                new_val = str(float(val.strip("%"))/100)
                if len(new_val.split(".")[1]) <= 2:
                    new_values.append(new_val)
                else:
                    new_values.append(str(round(100 * float(new_val))/100))
            return ",".join(new_values) 
            
    class Varscan_DepthTag():
        def __init__(self):
            self.metaheader = '##FORMAT=<ID=JQ_DP_VS,Number=1,Type=Float,Description="Jacquard depth for VarScan (based on DP),Source="Jacquard",Version={0}>\n'.format(jacquard_utils.__version__)
    
        def format(self, alt, filter, info_string, format_dict, count):
            if "DP" in format_dict.keys():
                format_dict["JQ_DP_VS"] = format_dict["DP"]
    
            return format_dict
        
    class Varscan_SomaticTag():
        def __init__(self):
            self.metaheader = '##FORMAT=<ID={0}VS,Number=1,Type=Integer,Description="Jacquard somatic status for VarScan: 0=non-somatic,1= somatic (based on SOMATIC info tag and if sample is TUMOR),Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_filter_tag, jacquard_utils.__version__)
    
        def format(self, alt, filter, info_string, format_dict, count):
            info_array = info_string.split(";")
            varscan_tag = jacquard_utils.jq_filter_tag + "VS"
    
            if "SS=2" in info_array:
                format_dict[varscan_tag] = self.somatic_status(count)
            else:
                format_dict[varscan_tag] = "0"
                
            return format_dict
    #  
        def somatic_status(self, count):
            if count == 0: #it's NORMAL
                return "0"
            else: #it's TUMOR
                return "1"

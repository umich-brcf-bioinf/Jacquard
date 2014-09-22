import jacquard.jacquard_utils as jacquard_utils
# import jacquard_utils


class Mutect():
    def __init__(self):
        self.name = "MuTect"
        
    def validate_input_file(self, input_file):
        valid = 0
        for line in input_file:
            if line.startswith("##MuTect"):
                valid = 1
            elif line.startswith("##"):
                continue
            else:
                break
        return (self.name, valid)

class AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}MT,Number=A,Type=Float,Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA),Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_af_tag, jacquard_utils.__version__)

    def format(self, alt, filter, info, format_dict, count):
        if "FA" in format_dict.keys():
            freq = format_dict["FA"].split(",")
            format_dict["JQ_AF_MT"] = self.roundTwoDigits(freq)
            
        return format_dict

    def roundTwoDigits(self, value): 
        new_values = []
        for val in value:
            if len(val.split(".")[1]) <= 2:
                new_values.append(val)
            else:
                new_values.append(str(round(100 * float(val))/100))
        return ",".join(new_values)
        
class DepthTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}MT,Number=1,Type=Float,Description="Jacquard depth for MuTect (based on DP),Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_dp_tag, jacquard_utils.__version__)

    def format(self, alt, filter, info, format_dict, count):
        if "DP" in format_dict.keys():
            format_dict["JQ_DP_MT"] = format_dict["DP"]
            
        return format_dict
    
class SomaticTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag),Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_somatic_tag, jacquard_utils.__version__)

    def format(self, alt, filter, info, format_dict, count):
        mutect_tag = jacquard_utils.jq_somatic_tag + "MT"
        if "SS" in format_dict.keys():
            format_dict[mutect_tag] = self.somatic_status(format_dict["SS"])
        else:
            format_dict[mutect_tag] = "0"
        return format_dict

    def somatic_status(self, ss_value):
        if ss_value == "2":
            return "1"
        else:
            return "0"

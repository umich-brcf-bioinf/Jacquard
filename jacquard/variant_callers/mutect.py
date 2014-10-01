import jacquard.jacquard_utils as jacquard_utils

class _AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}MT,Number=A,Type=Float,Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA),Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_af_tag, jacquard_utils.__version__)

    def format(self, vcfRecord):
        if "FA" in vcfRecord.format_set:
            sample_values = {}
            for key in vcfRecord.sample_dict.keys():
                freq = vcfRecord.sample_dict[key]["FA"].split(",")
                sample_values[key] = self._roundTwoDigits(freq)
            vcfRecord.insert_format_field("JQ_AF_MT",sample_values)

    def _roundTwoDigits(self, value): 
        new_values = []
        for val in value:
            if len(val.split(".")[1]) <= 2:
                new_values.append(val)
            else:
                new_values.append(str(round(100 * float(val))/100))
        return ",".join(new_values)
        
class _DepthTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}MT,Number=1,Type=Float,Description="Jacquard depth for MuTect (based on DP),Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_dp_tag, jacquard_utils.__version__)

    def format(self, vcfRecord):
        if "DP" in vcfRecord.format_set:
            sample_values = {}
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = vcfRecord.sample_dict[key]["DP"]
            vcfRecord.insert_format_field("JQ_DP_MT",sample_values)
    
class _SomaticTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag),Source="Jacquard",Version={1}>\n'.format(jacquard_utils.jq_somatic_tag, jacquard_utils.__version__)
        self.good = True
        
    def format(self, vcfRecord):
        mutect_tag = jacquard_utils.jq_somatic_tag + "MT"
        sample_values = {}
        if "SS" in vcfRecord.format_set:
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = self._somatic_status(vcfRecord.sample_dict[key]["SS"])
        else:
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = "0"
        vcfRecord.insert_format_field(mutect_tag,sample_values)  

    def _somatic_status(self, ss_value):
        if ss_value == "2":
            return "1"
        else:
            return "0"

class Mutect():
    def __init__(self):
        self.name = "MuTect"
        self.tags = [_AlleleFreqTag(),_DepthTag(),_SomaticTag()]
        
    def validate_input_file(self, header):
        valid = 0
        if "##MuTect" in header:
            valid = 1
        return (valid)
                
    def add_tags(self,vcfRecord):
        for tag in self.tags:
            tag.format(vcfRecord)
        return vcfRecord.asText()
    
    def update_metaheader(self,metaheader):
        for tag in self.tags:
            metaheader += tag.metaheader
        return metaheader    
    

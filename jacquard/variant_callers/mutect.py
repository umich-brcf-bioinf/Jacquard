import jacquard.utils as utils

class _AlleleFreqTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}MT,Number=A,Type=Float,Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA)",Source="Jacquard",Version={1}>'.format(utils.jq_af_tag, utils.__version__)

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
        self.metaheader = '##FORMAT=<ID={0}MT,Number=1,Type=Float,Description="Jacquard depth for MuTect (based on DP)",Source="Jacquard",Version={1}>'.format(utils.jq_dp_tag, utils.__version__)

    def format(self, vcfRecord):
        if "DP" in vcfRecord.format_set:
            sample_values = {}
            for key in vcfRecord.sample_dict.keys():
                sample_values[key] = vcfRecord.sample_dict[key]["DP"]
            vcfRecord.insert_format_field("JQ_DP_MT",sample_values)
    
class _SomaticTag():
    def __init__(self):
        self.metaheader = '##FORMAT=<ID={0}MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1=somatic (based on SS FORMAT tag)",Source="Jacquard",Version={1}>'.format(utils.jq_somatic_tag, utils.__version__)
        self.good = True
        
    def format(self, vcfRecord):
        mutect_tag = utils.jq_somatic_tag + "MT"
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
        
    def normalize(self, file_writer, file_readers):
        if len(file_readers) != 1:
                raise utils.JQException("ERROR: MuTect directories should have exactly one input file per patient, but found [{}].".format(len(file_readers)))

        file_writer.open()
        for file_reader in file_readers:
            file_reader.open()
            for line in file_reader.read_lines():
                file_writer.write(line)
            file_reader.close()
        file_writer.close()

    def get_new_metaheaders(self):
        return [tag.metaheader for tag in self.tags]

    def validate_input_file(self, meta_headers, column_header):
        valid = 0
        for line in meta_headers:
            if "##MuTect" in line:
                valid = 1
                break
        return (valid)
                
    def add_tags(self,vcfRecord):
        for tag in self.tags:
            tag.format(vcfRecord)
        return vcfRecord.asText()

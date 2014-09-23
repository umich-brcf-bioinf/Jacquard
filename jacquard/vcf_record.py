from collections import OrderedDict

class VcfRecord(object):

    def __init__(self, vcf_line):
        vcf_fields = vcf_line.rstrip("\n").split("\t")
        self.chrom,self.pos,self.id,self.ref,self.alt,self.qual,self.filter,self.info,self.format = vcf_fields[0:9]
        self.samples = vcf_fields[9:]
        tags = self.format.split(":")
        self.format_set = set(tags)
        self.sample_dict = {}
        for i,sample in enumerate(self.samples):
            values = sample.split(":")
            self.sample_dict[i] = OrderedDict(zip(tags,values))
                        
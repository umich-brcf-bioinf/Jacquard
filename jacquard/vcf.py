from __future__ import print_function
from collections import OrderedDict
import os


import jacquard_utils as jacquard_utils
from jacquard_utils import log

#TODO cgates: add context management to open/close 
class VcfReader(object):
    def __init__(self, file_path, get_caller):
        self.file_path = file_path
        self.name = os.path.basename(file_path)
        self.file_reader = None
        (self.column_header, self._metaheaders) = self._read_headers()
        self.caller = self._initialize_caller(get_caller)

    @property
    def metaheaders(self):
        return list(self._metaheaders)

    def _initialize_caller(self, get_caller):
        try:
            caller = get_caller(self.metaheaders, self.column_header, self.name)
            log("DEBUG: VCF [{}] recognized by caller [{}]",
                 self.name, caller.name)
            return caller
        except jacquard_utils.JQException as ex:
            log("ERROR: Problem parsing [{}]:{}", self.name, ex)
            raise ex


    def _read_headers(self):
        metaheaders = []
        with open(self.file_path, "r") as vcf:
            for line in vcf:
                if line.startswith("##"):
                    metaheaders.append(line.rstrip("\n"))
                elif line.startswith("#"):
                    column_header = line.rstrip("\n")
                else:
                    break
        return column_header, metaheaders

    def vcf_records(self):
        for line in self.file_reader:
            if line.startswith("#"):
                continue
            yield VcfRecord(line)

    def open(self):
        self.file_reader = open(self.file_path, 'r')

    def close(self):
        self.file_reader.close()


class VcfRecord(object):
    def __init__(self, vcf_line):
        vcf_fields = vcf_line.rstrip("\n").split("\t")
        self.chrom,self.pos,self.id,self.ref,self.alt,self.qual,self.filter,self.info,self.format = vcf_fields[0:9]
        self.samples = vcf_fields[9:]
        tags = self.format.split(":")
        self.format_set = tags
        self.sample_dict = {}
        for i,sample in enumerate(self.samples):
            values = sample.split(":")
            self.sample_dict[i] = OrderedDict(zip(tags,values))

    def asText(self):
        stringifier = [self.chrom,self.pos,self.id,self.ref,self.alt,self.qual,self.filter,self.info,":".join(self.format_set)]
        for key in self.sample_dict:
            stringifier.append(":".join(self.sample_dict[key].values()))
        return "\t".join(stringifier) + "\n"

    def insert_format_field(self, fieldname, field_dict):
        if fieldname in self.format_set:
            raise KeyError
        self.format_set.append(fieldname)
        if (field_dict.keys() != self.sample_dict.keys()):
            raise KeyError()
        for key in self.sample_dict.keys():
            self.sample_dict[key][fieldname] = str(field_dict[key])


#TODO cgates: add context management to open/close 
class VcfWriter(object):
    def __init__(self, output_filepath):
        self.output_filepath = output_filepath
        self._file_writer = None

    def open(self):
        self._file_writer = open(self.output_filepath, "w")

    def write(self, text):
        return self._file_writer.write(text)

    def close(self):
        self._file_writer.close()


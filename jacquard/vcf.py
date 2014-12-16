# pylint: disable=C0111
from __future__ import print_function
from collections import OrderedDict
import os
import sys

import utils

class RecognizedVcfReader(object):
    def __init__(self,vcf_reader,caller):
        self._vcf_reader = vcf_reader
        self.caller = caller

    @property
    def column_header(self):
        return self._vcf_reader.column_header

    def close(self):
        return self._vcf_reader.close() 

    @property
    def file_name(self):
        return self._vcf_reader.file_name

    @property
    def input_filepath(self):
        return self._vcf_reader.input_filepath

    @property
    def metaheaders(self):
        return self._vcf_reader.metaheaders

    def open(self):
        return self._vcf_reader.open()

    def vcf_records(self):
        return self._vcf_reader.vcf_records()


#TODO cgates: add context management to open/close
class VcfReader(object):
    def __init__(self, file_reader):
        self.input_filepath = file_reader.input_filepath
        self.file_name = file_reader.file_name
        self._file_reader = file_reader
        (self.column_header, self._metaheaders) = self._read_headers()

    @property
    def metaheaders(self):
        return list(self._metaheaders)

    #TODO: VcfReader shouldn't do open and close unless user tells it to
    def _read_headers(self):
        metaheaders = []
        column_header = None
        try:
            self._file_reader.open()
            for line in self._file_reader.read_lines():
                if line.startswith("##"):
                    metaheaders.append(line.rstrip())
                elif line.startswith("#"):
                    column_header = line.rstrip()
                else:
                    break
        finally:
            self._file_reader.close()

        if not (column_header and metaheaders):
            raise utils.JQException("ERROR: [{}] is not a valid vcf. Missing "
                                    "column header or metaheaders."\
                                    .format(self.file_name))

        return column_header, metaheaders

    def vcf_records(self):
        for line in self._file_reader.read_lines():
            if line.startswith("#"):
                continue
            yield VcfRecord(line)

    def open(self):
        self._file_reader.open()

    def close(self):
        self._file_reader.close()


class VcfRecord(object):
    def __init__(self, vcf_line):
        vcf_fields = vcf_line.rstrip().split("\t")
        self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, \
            self.filter, self.info, self.format = vcf_fields[0:9]
        self._key = (self._str_as_int(self.chrom), self.chrom, 
                     self._str_as_int(self.pos), self.ref, self.alt)
        self.samples = vcf_fields[9:]
        self.key = self.chrom+"_"+self.pos+"_"+self.ref+"_"+self.alt

        tags = self.format.split(":")
        self.format_set = tags

        self.sample_dict = {}
        for i, sample in enumerate(self.samples):
            values = sample.split(":")
            self.sample_dict[i] = OrderedDict(zip(tags, values))

    def get_info_dict(self):
        info_list = self.info.split(";")
        info_dict = {}

        for key_value in info_list:
            if "=" in key_value:
                key,value = key_value.split("=")
                info_dict[key] = value
            else:
                info_dict[key_value] = key_value

        return info_dict

    def asText(self):
        stringifier = [self.chrom, self.pos, self.id, self.ref, self.alt,
                   self.qual, self.filter, self.info,
                   ":".join(self.format_set)]

        for key in self.sample_dict:
            stringifier.append(":".join(self.sample_dict[key].values()))

        return "\t".join(stringifier) + "\n"

    def insert_format_field(self, fieldname, field_dict):
        if fieldname in self.format_set:
            raise KeyError
        self.format_set.append(fieldname)

        if field_dict.keys() != self.sample_dict.keys():
            raise KeyError()
        for key in self.sample_dict.keys():
            self.sample_dict[key][fieldname] = str(field_dict[key])

    def __eq__(self, other):
        return ("^".join([self.chrom,
                          self.pos,
                          self.ref,
                          self.alt]) == other)
        
    def __hash__(self):
        return hash("^".join([self.chrom,
                          self.pos,
                          self.ref,
                          self.alt]))
    
    def _str_as_int(self, string):
        if "chr" in string:
            string = string.replace("chr","")
        try:
            return int(string)
            
        except:
            return sys.maxint
    
    def __cmp__(self, other):
        return cmp(self._key, other._key)

#TODO cgates: add context management to open/close
class FileWriter(object):
    def __init__(self, output_filepath):
        self.output_filepath = output_filepath
        self._file_writer = None

    def open(self):
        self._file_writer = open(self.output_filepath, "w")

    def write(self, text):
        return self._file_writer.write(text)

    def close(self):
        self._file_writer.close()

    def __eq__(self, other):
        return (isinstance(other, self.__class__)
            and self.__dict__ == other.__dict__)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.output_filepath)

class FileReader(object):
    def __init__(self, input_filepath):
        self.input_filepath = input_filepath
        self.file_name = os.path.basename(input_filepath)
        self._file_reader = None

    def open(self):
        self._file_reader = open(self.input_filepath,"r")

    def read_lines(self):
        for line in self._file_reader:
            yield line

    def close(self):
        self._file_reader.close()

    def __eq__(self, other):
        return (isinstance(other, self.__class__)
            and self.__dict__ == other.__dict__)

    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(self.input_filepath)


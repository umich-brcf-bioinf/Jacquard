# pylint: disable=missing-docstring
from __future__ import print_function, absolute_import
from collections import OrderedDict, defaultdict
import jacquard.utils as utils
import os
import re
import sys


class RecognizedVcfReader(object):
    '''VcfReader with recognized caller'''
    def __init__(self, vcf_reader, caller):
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


#TODO: (cgates): add context management to open/close
class VcfReader(object):
    '''Wraps a file reader, providing VCF metaheaders and records'''

    def _init_format_metaheaders(self):
        format_ids = {}
        for metaheader in self.metaheaders:
            format_tag = re.match("^##FORMAT=.*?[<,]ID=([^,>]*)", metaheader)
            if format_tag:
                format_id = format_tag.group(1)
                format_ids[format_id] = metaheader.strip()
        return format_ids

    def __init__(self, file_reader):
        self._file_reader = file_reader
        (self.column_header, self.metaheaders) = self._init_headers()
        self._format_metaheaders = self._init_format_metaheaders()
        self.sample_names = self._init_sample_names()
        self.qualified_sample_names = self._create_qualified_sample_names()

    @property
    def file_name(self):
        return self._file_reader.file_name

    @property
    def input_filepath(self):
        return self._file_reader.input_filepath

    @property
    def format_metaheaders(self):
        return dict(self._format_metaheaders)

    def _init_sample_names(self):
        sample_names = []
        column_fields = self.column_header.split("\t")
        if column_fields > 8:
            sample_names = column_fields[9:]

        return sample_names

    def _create_qualified_sample_names(self):
        patient_prefix = self.file_name.split(".")[0]
        qualified_names = []
        for sample_name in self.sample_names:
            qualified_names.append("|".join([patient_prefix, sample_name]))
        return qualified_names

    def _init_headers(self):
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

        return column_header, tuple(metaheaders)

    def vcf_records(self, qualified=False):
        #pylint: disable=line-too-long
        sample_names = self.qualified_sample_names if qualified else self.sample_names

        for line in self._file_reader.read_lines():
            if line.startswith("#"):
                continue
            yield VcfRecord.parse_record(line, sample_names)

    def open(self):
        self._file_reader.open()

    def close(self):
        self._file_reader.close()

## pylint: disable=too-many-instance-attributes
# Alas, something must encapsulate the myriad VCF fields.
class VcfRecord(object):
    EMPTY_SET = set()

    @classmethod
    def parse_record(cls, vcf_line, sample_names):
        vcf_fields = vcf_line.rstrip("\n").split("\t")
        chrom, pos, rid, ref, alt, qual, rfilter, info \
                = vcf_fields[0:8]
        sample_fields = []
        sample_tag_values = {}
        if len(vcf_fields) > 9:
            rformat = vcf_fields[8]
            sample_fields = vcf_fields[9:]
            sample_tag_values = VcfRecord._sample_tag_values(sample_names,
                                                             rformat,
                                                             sample_fields)
        return VcfRecord(chrom, pos, ref, alt,
                         rid, qual, rfilter, info,
                         sample_tag_values)

    @classmethod
    def _sample_tag_values(cls, sample_names, rformat, sample_fields):
        sample_tag_values = OrderedDict()
        tag_names = VcfRecord._format_list(rformat)
        for i, sample_field in enumerate(sample_fields):
            tag_values = sample_field.split(":") if sample_field else "."
            sample_tag_values[sample_names[i]] = OrderedDict(zip(tag_names,
                                                                 tag_values))
        return sample_tag_values

    @classmethod
    def _format_list(cls, rformat):
        if rformat and rformat != ".":
            return rformat.split(":")
        else:
            return []

    @classmethod
    def _str_as_int(cls, string):
        if "chr" in string:
            string = string.replace("chr", "")
        try:
            return int(string)
        except ValueError:
            return sys.maxint

##TODO (cgates): adjust info field to be stored as dict instead of string
#TODO (cgates): adjust vcf names to not collide with reserved python words
## pylint: disable=too-many-arguments, invalid-name
# Alas, something must encapsulate the myriad VCF fields.
#  Note that some VCF field names collide with reserved python names
# (e.g. id, filter, format).
    def __init__(self, chrom, pos, ref, alt,
                 vcf_id=".", qual=".", vcf_filter=".", info=".",
                 sample_tag_values=None):
        self.chrom = chrom
        self.pos = pos
        self.id = vcf_id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = vcf_filter
        self.info = info
        self.info_dict = self.get_info_dict()

        if sample_tag_values is None:
            self.sample_tag_values = {}
        else:
            self.sample_tag_values = sample_tag_values
        self._key = (VcfRecord._str_as_int(self.chrom), self.chrom,
                     VcfRecord._str_as_int(self.pos), self.ref, self.alt)

    @property
    def format_tags(self):
        tags = VcfRecord.EMPTY_SET
        if self.sample_tag_values:
            first_sample = self.sample_tag_values.keys()[0]
            tags = set(self.sample_tag_values[first_sample].keys())
        return tags

    def get_info_dict(self):
        info_list = self.info.split(";")
        info_dict = {}
        for key_value in info_list:
            if "=" in key_value:
                key, value = key_value.split("=")
                info_dict[key] = value
            else:
                info_dict[key_value] = key_value

        return info_dict

    def add_info_field(self, field):
        if field in self.info_dict:
            msg = "New info field [{}] already exists.".format(field)
            raise KeyError(msg)

        if "=" in field:
            key, value = field.split("=")
            self.info_dict[key] = value
        else:
            self.info_dict[field] = field

        self._join_info_fields()

    def _join_info_fields(self):
        info_fields = []

        if len(self.info_dict) > 1:
            self.info_dict.pop(".", None)
        for field, value in self.info_dict.items():
            if field == value:
                info_fields.append(value)
            else:
                info_fields.append("=".join([field, value]))

        self.info = ";".join(info_fields)

    def get_empty_record(self):
        return VcfRecord(chrom=self.chrom,
                         pos=self.pos,
                         ref=self.ref,
                         alt=self.alt)

    def _format_field(self):
        format_field = "."
        if self.sample_tag_values:
            first_sample = self.sample_tag_values.keys()[0]
            tag_names = self.sample_tag_values[first_sample].keys()
            if tag_names:
                format_field = ":".join(tag_names)
        return format_field

    def _sample_field(self, sample):
        tag_values = self.sample_tag_values[sample].values()
        if tag_values:
            return ":".join(tag_values)
        else:
            return "."

    def asText(self):
        stringifier = [self.chrom, self.pos, self.id, self.ref, self.alt,
                       self.qual, self.filter, self.info,
                       self._format_field()]

        for sample in self.sample_tag_values:
            stringifier.append(self._sample_field(sample))

        return "\t".join(stringifier) + "\n"

    def _samples_match(self, new_sample_values):
        #pylint: disable=line-too-long
        return set(new_sample_values.keys()) == set(self.sample_tag_values.keys())

    def filter_sample_tag_values(self, desired_tags):
        new_sample_tag_values = defaultdict(dict)
        for sample, tags in self.sample_tag_values.items():
            for tag, value in tags.items():
                for desired_tag in desired_tags:
                    if re.match(desired_tag, tag):
                        new_sample_tag_values[sample][tag] = value

        self.sample_tag_values = new_sample_tag_values

    def add_sample_tag_value(self, tag_name, new_sample_values):
        if tag_name in self.format_tags:
            msg = "New format value [{}] already exists.".format(tag_name)
            raise KeyError(msg)

        if not self._samples_match(new_sample_values):
            raise KeyError("Sample name values must match "
                           "existing sample names")
        for sample in self.sample_tag_values.keys():
            value = str(new_sample_values[sample])
            self.sample_tag_values[sample][tag_name] = value

    def __eq__(self, other):
        return isinstance(other, VcfRecord) and self._key == other._key

    def __hash__(self):
        return hash(self._key)

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
        self._file_reader = open(self.input_filepath, "r")

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


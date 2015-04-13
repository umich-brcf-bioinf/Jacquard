"""Classes to parse, interpret, and manipulate VCF files and records."""
from __future__ import print_function, absolute_import, division

from collections import OrderedDict
import os
import re
import sys

import natsort

import jacquard.utils.utils as utils


#TODO: (cgates): add context management to open/close
class VcfReader(object):
    """Read only wrapper providing VCF metaheaders and records."""
    def __init__(self, file_reader):
        self._file_reader = file_reader
        (self.column_header, self.metaheaders) = self._init_headers()
        self.split_column_header = self.column_header.strip("#").split("\t")
        self.sample_names = self._init_sample_names()
        self.qualified_sample_names = self._create_qualified_sample_names()

    def __lt__(self, other):
        return self._file_reader < other._file_reader

    def _get_tag_metaheaders(self, regex_exp):
        tag_dict = {}
        for metaheader in self.metaheaders:
            tag = re.match(regex_exp, metaheader)
            if tag:
                tag_key = tag.group(1)
                tag_dict[tag_key] = metaheader.strip()

        return tag_dict

    @property
    def file_name(self):
        return self._file_reader.file_name

    @property
    def input_filepath(self):
        return self._file_reader.input_filepath

    @property
    def format_metaheaders(self):
        return dict(self._get_tag_metaheaders("^##FORMAT=.*?[<,]ID=([^,>]*)"))

    @property
    def info_metaheaders(self):
        return dict(self._get_tag_metaheaders("^##INFO=.*?[<,]ID=([^,>]*)"))

    @property
    def filter_metaheaders(self):
        return dict(self._get_tag_metaheaders("^##FILTER=.*?[<,]ID=([^,>]*)"))

    @property
    def contig_metaheaders(self):
        return dict(self._get_tag_metaheaders("^##contig=.*?[<,]ID=([^,>]*)"))

    def _init_sample_names(self):
        sample_names = []
        column_fields = self.column_header.split("\t")
        if len(column_fields) > 8:
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

    #TODO (cgates): qualified is used by ONE invocation in merge. Can we
    #somehow make merge do this instead of universally complicating the method?
    def vcf_records(self, qualified=False):
        """Generates parsed VcfRecord objects.

        Typically called in a for loop to process each vcf record in a
        VcfReader. VcfReader must be opened in advanced and closed when
        complete. Skips all headers.

        Args:
            qualified: When True, sample names are prefixed with file name

        Returns:
            Parsed VcfRecord

        Raises:
            StopIteration: when reader is exhausted.
            TypeError: if reader is closed.
        """
        if qualified:
            sample_names = self.qualified_sample_names
        else:
            sample_names = self.sample_names

        for line in self._file_reader.read_lines():
            if line.startswith("#"):
                continue
            yield VcfRecord.parse_record(line, sample_names)

    def open(self):
        self._file_reader.open()

    def close(self):
        self._file_reader.close()

class VcfRecord(object): #pylint: disable=too-many-instance-attributes
    """Represents an specific variant record.

    All functionality for parsing, interpreting, or manipulating a VCF record
    should be contained in this class. Mutator methods are implemented to add
    fields/tags to existing record, but replace/remove fields/tags.

    Attributes and args typically mirror the VCF fields; a few are renamed to
    avoid collisions with Python-reserved words. Except where noted, all
    attributes are stored as strings.

    Attributes:
            chrom:
            pos:
            ref:
            alt:
            vcf_id:
            qual:
            vcf_filter:
            info_dict: dict of info fields key-value pairs (flags stored as
                key-None)
            sample_tag_values: dict of samples where each value is a dict of
                tag-values. All keys and values are stored as strings. Example:
                {sampleA : {AF:0.3, DP:42}, sampleB : {AF:0.1, DP:25}}
                Note that dict is a natural pythonic representation of sample
                tag values, it's sometimes helpful to think of sample_tag_values
                as a table of samples (rows) x format_tags (columns).
            _key: tuple that defines record equality and ordering.
    """
    _EMPTY_SET = set()
    _FILTERS_TO_REPLACE = set(["", ".", "pass"])

    @classmethod
    def parse_record(cls, vcf_line, sample_names):
        """Alternative constructor that parses VcfRecord from VCF string.

        Aspire to parse/represent the data such that it could be reliably
        round-tripped. (This nicety means INFO fields and FORMAT tags should be
        treated as ordered to avoid shuffling.)

        Args:
            vcf_line: the VCF variant record as a string; tab separated fields,
                trailing newlines are ignored. Must have at least 8 fixed fields
                (through INFO)
            sample_names: a list of sample name strings; these should match
                the VCF header column
        Returns:
            A mutable VcfRecord.
        """
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
        """Creates a sample dict of tag-value dicts for a single variant record.

        Args:
            sample_names: list of sample name strings.
            rformat: record format string (from VCF record).
            sample_fields: list of strings where each string is the ';'
                seperated format values for an individual sample.

        Returns:
            An dict of samples, where each key is a sample and each value
                is an dict of format-values. See attribute below for example.
                Will return '.' if no values for sampe field.
        """
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
            return sys.maxsize
#TODO: (cgates): Could we make filter an OrderedSet
#TODO: (cgates) adjust info field to be stored as dict only instead of string
#pylint: disable=too-many-arguments
    def __init__(self, chrom, pos, ref, alt,
                 vcf_id=".", qual=".", vcf_filter=".", info=".",
                 sample_tag_values=None):
        """Builds a mutable VcfRecord from constituent VCF fields.

        If building from a VCF record line, see VcfRecord.parse_record.
        """
        self.chrom = chrom
        self.pos = pos
        self.vcf_id = vcf_id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = vcf_filter
        self.info = info
        self.info_dict = self._init_info_dict()

        if sample_tag_values is None:
            self.sample_tag_values = OrderedDict()
        else:
            self.sample_tag_values = sample_tag_values
        self._key = self._build_key()

    def _build_key(self):
        return (VcfRecord._str_as_int(self.chrom),
                self.chrom,
                VcfRecord._str_as_int(self.pos),
                self.ref,
                self.alt)

    @property
    def format_tags(self):
        """Returns set of format tags."""
        tags = VcfRecord._EMPTY_SET
        if self.sample_tag_values:
            first_sample = list(self.sample_tag_values.keys())[0]
            tags = set(self.sample_tag_values[first_sample].keys())
        return tags

    def _init_info_dict(self):
        info_dict = OrderedDict()
        if self.info and self.info != ".":
            info_list = self.info.split(";")
            for key_value in info_list:
                if "=" in key_value:
                    key, value = key_value.split("=")
                    info_dict[key] = value
                else:
                    info_dict[key_value] = key_value
        return info_dict

    def add_info_field(self, field):
        """Adds new info field (flag or key=value pair).

        Args:
            field: String flag (e.g. "SOMATIC") or key-value ("NEW_DP=42")

        Raises:
            KeyError: if info field already exists
        """
        if field in self.info_dict:
            msg = "New info field [{}] already exists.".format(field)
            raise KeyError(msg)

        if "=" in field:
            key, value = field.split("=")
            self.info_dict[key] = value
        else:
            self.info_dict[field] = field

        self._join_info_fields()

    #TODO:(cgates): Remove info; all external calls should reference info_dict
    def _join_info_fields(self):
        """Updates info attribute from info dict."""
        if self.info_dict:
            info_fields = []
            if len(self.info_dict) > 1:
                self.info_dict.pop(".", None)
            for field, value in self.info_dict.items():
                if field == value:
                    info_fields.append(value)
                else:
                    info_fields.append("=".join([field, value]))
            self.info = ";".join(info_fields)
        else:
            self.info = "."

    #TODO cgates: move this to merge
    def get_empty_record(self):
        return VcfRecord(chrom=self.chrom,
                         pos=self.pos,
                         ref=self.ref,
                         alt=self.alt)

    def _format_field(self):
        """Returns string representation of format field."""
        format_field = "."
        if self.sample_tag_values:
            first_sample = list(self.sample_tag_values.keys())[0]
            tag_names = self.sample_tag_values[first_sample].keys()
            if tag_names:
                format_field = ":".join(tag_names)
        return format_field

    def _sample_field(self, sample):
        """Returns string representation of sample-format values.

        Raises:
            KeyError: if requested sample is not defined.
        """
        tag_values = self.sample_tag_values[sample].values()
        if tag_values:
            return ":".join(tag_values)
        else:
            return "."

    def text(self):
        "Returns tab-delimited, newline terminated string of VcfRecord."
        stringifier = [self.chrom, self.pos, self.vcf_id, self.ref, self.alt,
                       self.qual, self.filter, self.info,
                       self._format_field()]

        for sample in self.sample_tag_values:
            stringifier.append(self._sample_field(sample))

        return "\t".join(stringifier) + "\n"

    def _samples_match(self, new_sample_values):
        return set(new_sample_values.keys()) == \
                set(self.sample_tag_values.keys())

    def add_sample_tag_value(self, tag_name, new_sample_values):
        """Appends a new format tag-value for all samples.

        Args:
            tag_name: string tag name; must not already exist
            new_sample

        Raises:
            KeyError: if tag_name to be added already exists
        """
        if tag_name in self.format_tags:
            msg = "New format value [{}] already exists.".format(tag_name)
            raise KeyError(msg)

        if not self._samples_match(new_sample_values):
            raise KeyError("Sample name values must match "
                           "existing sample names")
        for sample in self.sample_tag_values.keys():
            value = str(new_sample_values[sample])
            self.sample_tag_values[sample][tag_name] = value

    def add_or_replace_filter(self, new_filter):
        """Replaces null or blank filter or adds filter to existing list."""
        if self.filter.lower() in self._FILTERS_TO_REPLACE:
            self.filter = new_filter
        elif new_filter not in self.filter.split(";"):
            self.filter = ";".join([self.filter,
                                    new_filter])

    def __eq__(self, other):
        return isinstance(other, VcfRecord) and self._key == other._key

    def __hash__(self):
        return hash(self._key)

    def __lt__(self, other):
        return self._key < other._key


#TODO cgates: add context management to open/close
class FileWriter(object):
    """Trivial wrapper around os file to expedite testing."""
    def __init__(self, output_filepath):
        self.output_filepath = output_filepath
        self._file_writer = None

    @property
    def file_name(self):
        return os.path.basename(self.output_filepath)

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
    """Trivial wrapper around os file to expedite testing/natural sorting."""

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

    def __lt__(self, other):
        key = natsort.natsort_keygen()
        return key(self.file_name) < key(other.file_name)


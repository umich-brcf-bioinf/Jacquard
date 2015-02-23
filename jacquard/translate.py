from __future__ import absolute_import

import glob
import os

from jacquard import __version__
import jacquard.logger as logger
import jacquard.utils as utils
from jacquard.variant_callers import variant_caller_factory
from jacquard.vcf import FileReader, FileWriter


JQ_OUTPUT_SUFFIX = "translatedTags"


class _ExcludeMalformedRef(object):
    #pylint: disable=too-few-public-methods
    _VALID_REF = set(list("ACGTNacgtn"))
    _TAG_ID = "JQ_EXCLUDE_MALFORMED_REF"

    def __init__(self):
        self.metaheader = ('##FILTER=<'
                           'ID={},'
                           #pylint: disable=line-too-long
                           'Description="The format of the reference value for this variant record does not comply with VCF standard.",'
                           'Source="Jacquard",'
                           'Version="{}">').format(self._TAG_ID, __version__)

    def _is_valid_ref(self, record):
        return set(list(record.ref)).issubset(self._VALID_REF)

    def add_tag_values(self, record):
        if not self._is_valid_ref(record):
            record.add_or_replace_filter(self._TAG_ID)

class _ExcludeMalformedAlt(object):
    #pylint: disable=too-few-public-methods
    _VALID_ALT = set(list("*.ACGTNacgtn,"))
    _TAG_ID = "JQ_EXCLUDE_MALFORMED_ALT"

    def __init__(self):
        self.metaheader = ('##FILTER=<'
                           'ID={},'
                           #pylint: disable=line-too-long
                           'Description="The the format of the alternate allele value for this variant record does not comply with VCF standard.",'
                           'Source="Jacquard",'
                           'Version={}>').format(self._TAG_ID, __version__)

    def _is_valid_alt(self, record):
        valid_characters = set(list(record.alt)).issubset(self._VALID_ALT)
        invalid_characters = set(list("*.")).issubset(set(list(record.alt)))

        return valid_characters and not invalid_characters

    def add_tag_values(self, record):
        if not self._is_valid_alt(record):
            record.add_or_replace_filter(self._TAG_ID)

class _ExcludeMissingAlt(object):
    #pylint: disable=too-few-public-methods
    _MISSING_ALT = "."
    _TAG_ID = "JQ_EXCLUDE_MISSING_ALT"

    def __init__(self):
        self.metaheader = ('##FILTER='
                           '<ID={},'
                           #pylint: disable=line-too-long
                           'Description="The alternate allele is missing for this variant record.",'
                           'Source="Jacquard",'
                           'Version={}>').format(self._TAG_ID, __version__)

    def add_tag_values(self, record):
        if record.alt == self._MISSING_ALT:
            record.add_or_replace_filter(self._TAG_ID)

def _mangle_output_filename(input_file):
    basename, extension = os.path.splitext(os.path.basename(input_file))
    return ".".join([basename, JQ_OUTPUT_SUFFIX, extension.strip(".")])

def _write_headers(reader, new_tags, execution_context, file_writer):
    headers = reader.metaheaders
    headers.extend(execution_context)
    for tag in new_tags:
        headers.append(tag.metaheader)
    headers.append(reader.column_header)

    file_writer.write("\n".join(headers) + "\n")

def _build_file_readers(input_dir):
    in_files = glob.glob(os.path.join(input_dir, "*"))
    file_readers = []

    for in_file in in_files:
        file_readers.append(FileReader(in_file))

    return file_readers

def _translate_files(trans_vcf_reader,
                     new_tags,
                     execution_context,
                     file_writer):
    try:
        trans_vcf_reader.open()
        file_writer.open()

        _write_headers(trans_vcf_reader,
                       new_tags,
                       execution_context,
                       file_writer)

        records = []
        for record in trans_vcf_reader.vcf_records():
            for tag in new_tags:
                tag.add_tag_values(record)
            records.append(record)

        for record in sorted(records):
            file_writer.write(record.asText())

    finally:
        trans_vcf_reader.close()
        file_writer.close()

#TODO: Edit this later to be appropriate for translate.py
def add_subparser(subparser):
    #pylint: disable=line-too-long
    parser = subparser.add_parser("translate", help="Accepts a directory of VCf results and creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to Jacquard-tagged VCFs. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def execute(args, execution_context):
    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)

    file_readers = _build_file_readers(input_dir)

    trans_vcf_readers = variant_caller_factory.claim(file_readers)
    if not trans_vcf_readers:
        message = ("Specified input directory [{0}] contains no VCF files."
                   "Check parameters and try again.").format(input_dir)
        raise utils.JQException(message)

    logger.info("Processing [{}] VCF file(s) from [{}]",
                len(trans_vcf_readers),
                args.input)

    new_tags = [_ExcludeMalformedRef(),
                _ExcludeMalformedAlt(),
                _ExcludeMissingAlt()]

    for trans_vcf_reader in trans_vcf_readers:
        new_filename = _mangle_output_filename(trans_vcf_reader.file_name)
        output_filepath = os.path.join(output_dir, new_filename)
        file_writer = FileWriter(output_filepath)
        _translate_files(trans_vcf_reader,
                         new_tags,
                         execution_context,
                         file_writer,)

    logger.info("Wrote [{}] VCF file(s)", len(trans_vcf_readers))

def report_prediction(args):
    input_dir = os.path.abspath(args.input)
    file_readers = _build_file_readers(input_dir)
    output_file_names = set()

    for reader in file_readers:
        mangled_fname =_mangle_output_filename(reader.file_name)
        extension = os.path.splitext(os.path.basename(mangled_fname))[1]
        if extension == ".vcf":
            output_file_names.add(mangled_fname)

    return output_file_names

def get_required_input_output_types():
    return ("directory", "directory")


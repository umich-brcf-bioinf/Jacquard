from __future__ import absolute_import

import collections
import glob
import os

from jacquard import __version__
import jacquard.logger as logger
import jacquard.utils as utils
from jacquard.variant_callers import variant_caller_factory
from jacquard.vcf import FileReader, FileWriter


JQ_OUTPUT_SUFFIX = "translatedTags"

#TODO: Edit this later to be appropriate for translate.py
def add_subparser(subparser):
    #pylint: disable=line-too-long
    parser = subparser.add_parser("translate", help="Accepts a directory of VCf results and creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to Jacquard-tagged VCFs. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def _mangle_output_filenames(input_file):
    basename, extension = os.path.splitext(os.path.basename(input_file))
    return ".".join([basename, JQ_OUTPUT_SUFFIX, extension.strip(".")])

def _get_output_filenames(input_files):
    output_files = set()
    for input_file in input_files:
        output_files.add(_mangle_output_filenames(input_file))

    return output_files

def _predict_output(args):
    input_dir = os.path.abspath(args.input)
    input_files = sorted(glob.glob(os.path.join(input_dir, "*.vcf")))
    desired_output_files = _get_output_filenames(input_files)

    return desired_output_files

def report_prediction(args):
    return _predict_output(args)

def get_required_input_output_types():
    return ("directory", "directory")

def _comma_separated(line):
    split_line = line.split(",")
    for char in split_line:
        if char == "":
            return False
    return True

def _check_records(reader):
    correct_ref = "ACGTNacgtn"
    correct_alt = "*.ACGTNacgtn,"
    anomalous_set = set()
    anomalous_records = []

    malformed_ref = 0
    malformed_alt = 0
    missing_alt = 0

    for vcf_record in reader.vcf_records():
        anomalous_flags = []
        #TODO: (cgates): Suspect there may be a clearer way to do this?
        if _comma_separated(vcf_record.ref) and not\
                    all(x in correct_ref for x in list(vcf_record.ref)):
            anomalous_flags.append("JQ_MALFORMED_REF")
            malformed_ref += 1
        if _comma_separated(vcf_record.alt) and not\
                    all(x in correct_alt for x in list(vcf_record.alt)):
            anomalous_flags.append("JQ_MALFORMED_ALT")
            malformed_alt += 1
        if vcf_record.alt == "*" or vcf_record.alt == ".":
            anomalous_flags.append("JQ_MISSING_ALT")
            missing_alt += 1
        if len(anomalous_flags) > 0:
            anomalous_flags = ["JQ_EXCLUDE"] + anomalous_flags

        for flag in anomalous_flags:
            anomalous_set.add(flag)
        anomalous_records.append(anomalous_flags)

    if malformed_ref:
        logger.debug("{}|Added filter flag [JQ_MALFORMED_REF] to [{}] variant "
                     "records.", reader.file_name, malformed_ref)
    if malformed_alt:
        logger.debug("{}|Added filter flag [JQ_MALFORMED_ALT] to [{}] variant "
                     "records.", reader.file_name, malformed_alt)
    if missing_alt:
        logger.debug("{}|Added filter flag [JQ_MISSING_ALT] to [{}] variant "
                     "records.", reader.file_name, missing_alt)

    return anomalous_set, anomalous_records

def _write_headers(file_writer, reader, execution_context):
    headers = reader.metaheaders
    headers.extend(execution_context)
    headers.append("##jacquard.tag.caller={0}".format(reader.caller_name))
    headers.append(reader.column_header)

    print headers
    file_writer.write("\n".join(headers) + "\n")

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
        if not self._is_valid_alt(record) :
            record.add_or_replace_filter(self._TAG_ID)

class _ExcludeMissingAlt(object):
    #pylint: disable=too-few-public-methods
    _VALID_ALT = set(list("*."))
    _TAG_ID = "JQ_EXCLUDE_MISSING_ALT"

    def __init__(self):
        self.metaheader = ('##FILTER='
                           '<ID={},'
                           #pylint: disable=line-too-long
                           'Description="The alternate allele is missing for this variant record.",'
                           'Source="Jacquard",'
                           'Version={}>').format(self._TAG_ID, __version__)

    def _is_valid_alt(self, record):
        valid_length = len(record.alt) == 1
        return set(list(record.alt)).issubset(self._VALID_ALT) and valid_length

    def add_tag_values(self, record):
        if self._is_valid_alt(record):
            record.add_or_replace_filter(self._TAG_ID)


def _build_file_readers(input_dir):
    in_files = glob.glob(os.path.join(input_dir, "*.vcf"))
    file_readers = []

    for in_file in in_files:
        file_readers.append(FileReader(in_file))

    return file_readers

#TODO: execution_context is undergoing construction.
def _translate_files(trans_vcf_reader, file_writer, execution_context):
    trans_vcf_reader.open()
    file_writer.open()

    trans_vcf_reader.add_tag_class([_ExcludeMalformedRef(),
                                    _ExcludeMalformedAlt(),
                                    _ExcludeMissingAlt()])

    _write_headers(file_writer, trans_vcf_reader, execution_context)
    for record in trans_vcf_reader.vcf_records():
        file_writer.write(record.asText())

    trans_vcf_reader.close()
    file_writer.close()

# def _translate_files(trans_vcf_readers,
#                      output_dir,
#                      execution_context):
#     total_filtered_records = 0
#     callers = collections.defaultdict(int)
# 
#     for trans_vcf_reader in trans_vcf_readers:
#         new_filename = _mangle_output_filenames(trans_vcf_reader.file_name)
#         output_filepath = os.path.join(output_dir, new_filename)
#         file_writer = FileWriter(output_filepath)
# 
#         trans_vcf_reader.open()
#         anomalous_set, anomalous_records = _check_records(trans_vcf_reader)
#         callers[trans_vcf_reader.caller_name] += len(anomalous_records)
#         trans_vcf_reader.close()
# 
#         trans_vcf_reader.open()
#         file_writer.open()
#         _write_headers(file_writer,
#                        trans_vcf_reader,
#                        execution_context,
#                        anomalous_set)
#         _write_records(file_writer,
#                        trans_vcf_reader.vcf_records(),
#                        anomalous_records)
#         trans_vcf_reader.close()
#         file_writer.close()
# 
#     for caller in callers:
#         total_filtered_records += callers[caller]
#         if callers[caller]:
#             logger.debug("Added a filter flag to [{}] problematic {} variant "
#                          "records.", callers[caller], caller)
# 
#     if total_filtered_records:
#         logger.warning("A total of [{}] problematic variant records failed "
#                        "Jacquard's filters. See output and log for details.",
#                        total_filtered_records)

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

    for trans_vcf_reader in trans_vcf_readers:
        new_filename = _mangle_output_filenames(trans_vcf_reader.file_name)
        output_filepath = os.path.join(output_dir, new_filename)
        file_writer = FileWriter(output_filepath)
        _translate_files(trans_vcf_reader, file_writer, execution_context)

    logger.info("Wrote [{}] VCF file(s)", len(trans_vcf_readers))

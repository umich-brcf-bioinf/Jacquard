"""Translates a set of VCF files by adding standardized tags.

Reads incoming VCF files determining appropriate "origin caller" (e.g. MuTect);
emits a new translated file. The translated file is similar to the input
file, with these exceptions:
    * translate will add a filter flag anomalous VCF records, i.e. records
        that don't conform to the standard; for example both Strelka and
        Varscan emit VCF records with invalid ALT values.
    * translate will add new Jacquard-standard FORMAT tags that augment the
        caller specific tags (e.g. a Varscan FREQ tag would generates a new
        JQ_VS_AF tag).

There will typically be a translated VCF file for each input VCF file.
Unrecognized VCFs are not copied to output.
The origin caller of a VCFs is recognized in part by the content of metaheaders,
so it's imperative that VCF metaheaders be present and accurate.
"""
from __future__ import print_function, absolute_import, division

import argparse
from collections import defaultdict
import glob
import os

import jacquard.utils.logger as logger
import jacquard.utils.utils as utils
from jacquard.variant_caller_transforms import variant_caller_factory
from jacquard.utils.vcf import FileReader, FileWriter


_FILE_OUTPUT_SUFFIX = "translatedTags"


class _ExcludeMalformedRef(object):
    #pylint: disable=too-few-public-methods
    _VALID_REF = set(list("ACGTNacgtn"))
    _TAG_ID = "JQ_EXCLUDE_MALFORMED_REF"

    def __init__(self):
        self.metaheader = ('##FILTER=<'
                           'ID={},'
                           #pylint: disable=line-too-long
                           'Description="The format of the reference value for this variant record does not comply with VCF standard.">')\
                           .format(self._TAG_ID)

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
                           'Description="The the format of the alternate allele value for this variant record does not comply with VCF standard.">')\
                           .format(self._TAG_ID)

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
                           'Description="The alternate allele is missing for this variant record.">')\
                           .format(self._TAG_ID)

    def add_tag_values(self, record):
        if record.alt == self._MISSING_ALT:
            record.add_or_replace_filter(self._TAG_ID)


def _build_file_readers(input_dir):
    in_files = glob.glob(os.path.join(input_dir, "*"))
    file_readers = []

    for in_file in in_files:
        file_readers.append(FileReader(in_file))

    return file_readers

def _mangle_output_filename(input_file):
    basename, extension = os.path.splitext(os.path.basename(input_file))
    return ".".join([basename, _FILE_OUTPUT_SUFFIX, extension.strip(".")])

def _check_snp_indel_pairings(altered_file_names, args):
    if not set([len(i) for i in altered_file_names.values()]) == set([1]):
        if not args.allow_inconsistent_sample_sets:
            error = 0
            for file_names in altered_file_names.values():
                if len(file_names) % 2 != 0:
                    message = ("File {} was missing a "
                               "corresponding snp/indel file.")
                    logger.error(message, file_names)
                    error = 1
            if error:
                message = ("Some VCFs were missing either a snp/snvs "
                            "or an indel/indels file. Review "
                            "inputs/command options to align file "
                            "pairings or use the flag "
                            "--allow_inconsistent_sample_sets.")
                raise utils.UsageError(message)

def _validate_file_pairings(args, claimed_vcf_readers):
    input_path = args.input
    if not os.path.isdir(input_path):
        input_path = os.path.dirname(input_path)

    input_vcfs = sorted(glob.glob(os.path.join(input_path, "*.vcf")))
    altered_file_names = defaultdict(list)
    for vcf_reader in claimed_vcf_readers:
        for input_vcf in input_vcfs:
            basename = os.path.basename(input_vcf)
            word_list = vcf_reader.expected_file_format()
            file_names = [i for i in basename.split(".") if i not in word_list]
            joined_file_names = ".".join(file_names)
            if len(joined_file_names) != len(basename):
                altered_file_names[joined_file_names].append(basename)

    _check_snp_indel_pairings(altered_file_names, args)

def validate_args(args):
    unclaimed_readers, trans_vcf_readers = _claim_readers(args)
    if unclaimed_readers and not args.force:
        raise utils.UsageError(_build_validation_message(unclaimed_readers))
    elif not trans_vcf_readers:
        raise utils.UsageError(("no vcfs in input dir "
                                "[{}] can be translated.").format(args.input))
    _validate_file_pairings(args, trans_vcf_readers)

def _build_validation_message(unclaimed_readers):
    total_unclaimed = len(unclaimed_readers)
    if total_unclaimed == 1:
        unclaimed_details = "file [{}]".format(unclaimed_readers[0].file_name)
    else:
        cutoff = 5
        file_names = [reader.file_name for reader in unclaimed_readers]
        file_names = file_names[0:min(cutoff, total_unclaimed)]
        unclaimed_list = ", ".join(file_names)
        if total_unclaimed > cutoff:
            omitted = total_unclaimed - cutoff
            unclaimed_list += ", ...({} file(s) omitted)".format(omitted)
        unclaimed_details = "files [{}]".format(unclaimed_list)

    return ("{} input {} cannot be "
            "translated; review input dir and try again, or "
            "use '--force' to ignore these files.").format(total_unclaimed,
                                                           unclaimed_details)

def _claim_readers(args):
    input_dir = os.path.abspath(args.input)
    file_readers = _build_file_readers(input_dir)
    factory = variant_caller_factory.VariantCallerFactory(args)

    return factory.claim(file_readers)

def _log_unclaimed_readers(unclaimed_readers):
    unclaimed_log_message = "The input file [{}] will not be translated"
    for reader in sorted(unclaimed_readers):
        msg = unclaimed_log_message.format(reader.file_name)
        logger.warning(msg)

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

        for record in trans_vcf_reader.vcf_records():
            for tag in new_tags:
                tag.add_tag_values(record)
            file_writer.write(record.text())

    finally:
        trans_vcf_reader.close()
        file_writer.close()

def _write_headers(reader, new_tags, execution_context, file_writer):
    headers = reader.metaheaders
    headers.extend(execution_context)
    for tag in new_tags:
        headers.append(tag.metaheader)

    sorted_headers = utils.sort_metaheaders(headers)
    sorted_headers.append(reader.column_header)

    file_writer.write("\n".join(sorted_headers) + "\n")

def get_required_input_output_types():
    return ("directory", "directory")

def report_prediction(args):
    input_dir = os.path.abspath(args.input)
    file_readers = _build_file_readers(input_dir)
    output_file_names = set()

    for reader in file_readers:
        mangled_fname = _mangle_output_filename(reader.file_name)
        extension = os.path.splitext(os.path.basename(mangled_fname))[1]
        if extension == ".vcf":
            output_file_names.add(mangled_fname)

    return output_file_names

def add_subparser(subparser):
    #pylint: disable=line-too-long
    parser = subparser.add_parser("translate", formatter_class=argparse.RawTextHelpFormatter, help="Accepts a directory of VCF results (and VarScan high confidence files). Creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser.add_argument("input", help="Path to directory containing VCFs (and VarScan high confidence files). Other file types ignored")
    parser.add_argument("output", help="Path to Jacquard-tagged VCFs. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")
    parser.add_argument("--varscan_hc_filter_file_regex",
                        help=("Regex pattern that identifies optional VarScan high-confidence filter files.\n"
                              "The VCF, high-confidence file pairs should share the same prefix.\n"
                              "For example, given patientA.snp.vcf, patientA.indel.vcf, patientA.snp.fpfilter.pass, patientA.indel.fpfilter.pass,\n"
                              "you could enable this option as varscan_hc_filter_file_regex='.fpfilter.pass$'"))
    parser.add_argument("--allow_inconsistent_sample_sets", help="Allow inconsistent sample sets to be used")
    parser.add_argument("--log_file", help="Path to log file destination. Defaults to current working directory if not specified.")

#TODO (cgates): This module is both a command and also manipulates VcfRecords
# like a caller. This is the only body of code that does both these things.
# Does this bother anyone else?
def execute(args, execution_context):
    validate_args(args)

    output_dir = os.path.abspath(args.output)
    unclaimed_readers, trans_vcf_readers = _claim_readers(args)

    _log_unclaimed_readers(unclaimed_readers)

    logger.info("Processing [{}] VCF file(s) from [{}]",
                len(trans_vcf_readers),
                args.input)

    new_tags = [_ExcludeMalformedRef(),
                _ExcludeMalformedAlt(),
                _ExcludeMissingAlt()]

    for i, trans_vcf_reader in enumerate(trans_vcf_readers):
        logger.info("Translating file {}/{} [{}]",
                    i + 1,
                    len(trans_vcf_readers),
                    trans_vcf_reader.file_name)
        new_filename = _mangle_output_filename(trans_vcf_reader.file_name)
        output_filepath = os.path.join(output_dir, new_filename)
        file_writer = FileWriter(output_filepath)
        _translate_files(trans_vcf_reader,
                         new_tags,
                         execution_context,
                         file_writer,)

    logger.info("Wrote [{}] VCF file(s)", len(trans_vcf_readers))

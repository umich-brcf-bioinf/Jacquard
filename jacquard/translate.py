from __future__ import absolute_import
from jacquard import __version__
from jacquard.variant_callers import variant_caller_factory
from jacquard.vcf import FileReader, FileWriter
import glob
import jacquard.logger as logger
import jacquard.utils as utils
import os




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

        for record in trans_vcf_reader.vcf_records():
            for tag in new_tags:
                tag.add_tag_values(record)
            file_writer.write(record.asText())

    finally:
        trans_vcf_reader.close()
        file_writer.close()

def add_subparser(subparser):
    #pylint: disable=line-too-long
    parser = subparser.add_parser("translate", help="Accepts a directory of VCF results (and VarScan high confidence files). Creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser.add_argument("input", help="Path to directory containing VCFs (and VarScan high confidence files). Other file types ignored")
    parser.add_argument("output", help="Path to Jacquard-tagged VCFs. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")


def _log_unclaimed_readers(unclaimed_readers):
    unclaimed_log_messgae = "The input file [{}] will not be translated"
    for reader in unclaimed_readers:
        msg = unclaimed_log_messgae.format(reader.file_name)
        logger.warning(msg)

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
        mangled_fname = _mangle_output_filename(reader.file_name)
        extension = os.path.splitext(os.path.basename(mangled_fname))[1]
        if extension == ".vcf":
            output_file_names.add(mangled_fname)

    return output_file_names

def get_required_input_output_types():
    return ("directory", "directory")

def _claim_readers(args):
    input_dir = os.path.abspath(args.input)
    file_readers = _build_file_readers(input_dir)
    return variant_caller_factory.claim(file_readers)

def validate_args(args):
    unclaimed_readers, trans_vcf_readers = _claim_readers(args)
    if unclaimed_readers and not args.force:
        raise utils.UsageError(_build_validation_message(unclaimed_readers))
    elif not trans_vcf_readers:
        raise utils.UsageError(("no vcfs in input dir "
                                "[{}] can be translated.").format(args.input))

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


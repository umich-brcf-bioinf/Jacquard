from __future__ import print_function, absolute_import
import collections
import glob
import os
import shutil

from jacquard.variant_callers import variant_caller_factory
import jacquard.utils as utils
import jacquard.vcf as vcf
import jacquard.logger as logger

#pylint: disable=line-too-long
def add_subparser(subparser):
    parser = subparser.add_parser("tag", help="Accepts a directory of VCf results and creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to Jacquard-tagged VCFs. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def _comma_separated(line):
    split_line = line.split(",")
    for char in split_line:
        if char == "":
            return False
    return True

#TODO: (cgates): This works, but would rather see this as a collection of validator methods
def _check_records(reader):
    correct_ref = "ACGTNacgtn"
    correct_alt = "*.ACGTNacgtn"
    anomalous_set = set()
    anomalous_records = []

    malformed_ref = 0
    malformed_alt = 0
    missing_alt = 0

    for vcf_record in reader.vcf_records():
        anomalous_flags = []
        if _comma_separated(vcf_record.ref) and not all(x in correct_ref for x in list(vcf_record.ref)):
            anomalous_flags.append("JQ_MALFORMED_REF")
            malformed_ref += 1
        if _comma_separated(vcf_record.alt) and not all(x in correct_alt for x in list(vcf_record.alt)):
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

def _write_records(reader, writer, anomalous_records):
    for i, vcf_record in enumerate(reader.vcf_records()):
        anomalous = ";".join(anomalous_records[i])
        if anomalous:
            if vcf_record.filter == ".":
                vcf_record.filter = anomalous
            else:
                vcf_record.filter = ";".join([vcf_record.filter, anomalous])

        writer.write(reader.caller.add_tags(vcf_record))

def _modify_metaheaders(reader, writer, execution_context, anomalous_set):
    new_headers = list(reader.metaheaders)
    new_headers.extend(execution_context)
    new_headers.append("##jacquard.tag.caller={0}".format(reader.caller.name))
    new_headers.extend(reader.caller.get_new_metaheaders())

    if len(anomalous_set) > 0:
        new_headers.append('##FILTER=<ID=JQ_EXCLUDE,Description="This ' \
                           'variant record is problematic and will be '\
                           'excluded from downstream Jacquard processing.",'
                           'Source="Jacquard",Version="">')
        if "JQ_MALFORMED_REF" in anomalous_set:
            new_headers.append('##FILTER=<ID=JQ_MALFORMED_REF,Description='\
                               '"The format of the reference value for this '\
                               'variant record does not comply with VCF '\
                               'standard.",Source="Jacquard",Version="">')
        if "JQ_MALFORMED_ALT" in anomalous_set:
            new_headers.append('##FILTER=<ID=JQ_MALFORMED_ALT,Description='\
                               '"The the format of the alternate allele value '\
                               'for this variant record does not comply with '\
                               'VCF standard.",Source="Jacquard",Version="">')
        if "JQ_MISSING_ALT" in anomalous_set:
            new_headers.append('##FILTER=<ID=JQ_MISSING_ALT,Description="The '\
                               'alternate allele is missing for this variant '\
                               'record.",Source="Jacquard",Version="">')
    new_headers.append(reader.column_header)

    writer.write("\n".join(new_headers) +"\n")

def tag_files(vcf_readers_to_writers, execution_context,
              get_caller=variant_caller_factory.get_caller):

    total_number_of_files = len(vcf_readers_to_writers)
    callers = collections.defaultdict(int)
    for count, item in enumerate(vcf_readers_to_writers.items()):
        reader, writer = item
        logger.info("Reading [{}] ({}/{})",
                    reader.input_filepath,
                    count + 1,
                    total_number_of_files)
        reader.open()

        anomalous_set, anomalous_records = _check_records(reader)
        caller = get_caller(reader.metaheaders,
                            reader.column_header,
                            reader.file_name)
        callers[caller.name] += len(anomalous_records)

        reader.close()

        reader.open()
        writer.open()
        _modify_metaheaders(reader, writer, execution_context, anomalous_set)
        _write_records(reader, writer, anomalous_records)
        reader.close()
        writer.close()

    total_filtered_records = 0

    for caller in callers:
        total_filtered_records += callers[caller]
        if callers[caller]:
            logger.debug("Added a filter flag to [{}] problematic {} variant "
                         "records.", callers[caller], caller)

    if total_filtered_records:
        logger.warning("A total of [{}] problematic variant records failed "
                       "Jacquard's filters. See output and log for details.",
                       total_filtered_records)

def _log_caller_info(vcf_readers):
    caller_count = collections.defaultdict(int)
    for vcf_reader in vcf_readers:
        caller_count[vcf_reader.caller.name] += 1
    for caller_name in sorted(caller_count):
        logger.info("Recognized [{}] {} file(s)",
                    caller_count[caller_name], caller_name)

def _build_vcf_readers(input_dir,
                       get_caller=variant_caller_factory.get_caller):
    in_files = sorted(glob.glob(os.path.join(input_dir, "*.vcf")))
    vcf_readers = []
    failures = 0
    for filename in in_files:
        file_path = os.path.join(input_dir, filename)
        try:
            vcf_reader = vcf.VcfReader(vcf.FileReader(file_path))
            caller = get_caller(vcf_reader.metaheaders,
                                vcf_reader.column_header,
                                vcf_reader.file_name)
            vcf_readers.append(vcf.RecognizedVcfReader(vcf_reader, caller))
        except utils.JQException:
            failures += 1
    if failures:
        raise utils.JQException(("[{}] VCF files could not be parsed."
                                 " Review logs for details, adjust input, "
                                 "and try again.").format(failures))
    _log_caller_info(vcf_readers)
    return vcf_readers


def _build_vcf_readers_to_writers(vcf_readers, output_dir):
    vcf_providers_to_writers = {}
    for reader in vcf_readers:
        basename, extension = os.path.splitext(reader.file_name)
        new_filename = basename + ".jacquardTags" + extension
        output_filepath = os.path.join(output_dir, new_filename)
        vcf_providers_to_writers[reader] = vcf.FileWriter(output_filepath)

    return vcf_providers_to_writers

def execute(args, execution_context):
    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    utils.validate_directories(input_dir, output_dir)

    vcf_readers = _build_vcf_readers(input_dir)
    if not vcf_readers:
        logger.error(("Specified input directory [{0}] contains no VCF files."
                      "Check parameters and try again."),
                     input_dir)
        #TODO cgates: move to jacquard.py
        shutil.rmtree(output_dir)
        exit(1)

    readers_to_writers = _build_vcf_readers_to_writers(vcf_readers, output_dir)
    logger.info("Processing [{}] VCF file(s) from [{}]",
                len(vcf_readers),
                input_dir)

    tag_files(readers_to_writers, execution_context)
    logger.info("Wrote [{}] VCF file(s) to [{}]",
                len(readers_to_writers),
                output_dir)

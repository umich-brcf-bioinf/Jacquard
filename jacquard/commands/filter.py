#pylint: disable=too-many-locals
from __future__ import print_function, absolute_import

from collections import defaultdict
import collections
import glob
import os
import re

import natsort

import jacquard.logger as logger
import jacquard.utils as utils
import jacquard.variant_caller_transforms.variant_caller_factory as variant_caller_factory
import jacquard.vcf as vcf


_FILE_OUTPUT_SUFFIX = "HCsomatic"
_JQ_SOMATIC_TAG = "HC_SOM"

#TODO: (cgates): Use tuples instead of concatenated strings
#TODO: (cgates): refactor this as a stats object that collects info in
# the main processing loop

def _find_passed_variants(record, positions, num_filtered_records):
    if "PASS" in record.filter:
        passed_key = "^".join([record.chrom,
                               record.pos,
                               record.ref])
        positions[passed_key] = 1
    else:
        num_filtered_records += 1

    return positions, num_filtered_records

def _find_som_variants(record, positions, num_filtered_records):
    sample_tag_values = record.sample_tag_values
    somatic = 0
    for sample in sample_tag_values:
        for tag in sample_tag_values[sample]:
            if re.search(_JQ_SOMATIC_TAG, tag):
                if sample_tag_values[sample][tag] == "1":
                    somatic = 1
                    somatic_key = "^".join([record.chrom,
                                            record.pos,
                                            record.ref])
                    positions[somatic_key] = 1
    if somatic == 0:
        num_filtered_records += 1

    return positions, num_filtered_records

def _handle_include_variants(vcf_reader, positions, include_variants):
    num_filtered_records = 0
    vcf_reader.open()

    for record in vcf_reader.tagged_vcf_records():
        if include_variants == "passed":
            (positions,
             num_filtered_records) = _find_passed_variants(record,
                                                           positions,
                                                           num_filtered_records)
        elif include_variants == "somatic":
            (positions,
             num_filtered_records) = _find_som_variants(record,
                                                        positions,
                                                        num_filtered_records)

    vcf_reader.close()

    return num_filtered_records, positions

def _handle_include_loci(vcf_reader, positions, include_loci):
    num_filtered_records = 0
    vcf_reader.open()

    for record in vcf_reader.tagged_vcf_records():
        if include_loci == "all_passed" or include_loci == "any_passed":
            (positions,
             num_filtered_records) = _find_passed_variants(record,
                                                           positions,
                                                           num_filtered_records)

        elif include_loci == "all_somatic" or include_loci == "any_somatic":
            (positions,
            num_filtered_records) = _find_som_variants(record,
                                                       positions,
                                                       num_filtered_records)

    vcf_reader.close()

    return num_filtered_records, positions

def _get_default_positions(vcf_reader, positions):
    num_filtered_records = 0
    vcf_reader.open()

    for record in vcf_reader.tagged_vcf_records():
        if "JQ_EXCLUDE" in record.filter:
            num_filtered_records += 1
        else:
            (positions,
            num_filtered_records) = _find_som_variants(record,
                                                       positions,
                                                       num_filtered_records)
    vcf_reader.close()

    return num_filtered_records, positions

def _find_positions(vcf_readers, include_variants, include_loci):
    positions = {}
    no_jq_tags = []
    filtered_per_caller = defaultdict(int)

    for i, vcf_reader in enumerate(vcf_readers):
        logger.info("Reading [{}] ({}/{})",
                    os.path.basename(vcf_reader.file_name),
                    i + 1,
                    len(vcf_readers))
        if include_variants:
            (num_filtered_records,
             positions) = _handle_include_variants(vcf_reader,
                                                   positions,
                                                   include_variants)
        elif include_loci:
            (num_filtered_records,
             positions) = _handle_include_loci(vcf_reader,
                                               positions,
                                               include_loci)
        else:
            num_filtered_records, positions = _get_default_positions(vcf_reader,
                                                                     positions)

        filtered_per_caller[vcf_reader.caller_name] += num_filtered_records

        if len(positions) == 0:
            no_jq_tags.append(vcf_reader.file_name)
            logger.warning("Input file [{}] has no high-confidence somatic "
                           "variants.", os.path.basename(vcf_reader.file_name))
    _log_positions(positions, filtered_per_caller, no_jq_tags)

    return positions

def _log_positions(somatic_positions, filtered_per_caller, no_jq_tags):
    if no_jq_tags:
        logger.warning("[{}] VCF file(s) had no high-confidence somatic "
                       "variants. See log for details.", len(no_jq_tags))

    for caller, num_filtered in list(filtered_per_caller.items()):
        if num_filtered > 0:
            logger.debug("Removed [{}] problematic {} variant records with "
                         "filter=JQ_EXCLUDE",
                         num_filtered,
                         caller)

    total_filtered_records = sum(filtered_per_caller.values())
    if total_filtered_records > 0:
        logger.warning("A total of [{}] problematic variant records failed "
                       "Jacquard's filters. See output and log for details.",
                       total_filtered_records)

    total_somatic_records = sum(somatic_positions.values())
    total_searched_records = total_somatic_records + total_filtered_records
    logger.info("Searched [{}] variant calls.", total_searched_records)

    logger.info("Found [{}] distinct high-confidence somatic positions.",
                len(somatic_positions.keys()))

def _add_headers(execution_context, somatic_positions):
    somatic_positions_header = "##jacquard.filterHCSomatic."\
                               "total_highConfidence_somatic_positions={0}\n"\
                               .format(len(somatic_positions.keys()))

    if not execution_context:
        execution_context = []
    execution_context.append(somatic_positions_header)

    return execution_context

def _write_positions(readers_to_writers,
                     positions,
                     execution_context,
                     include_loci):

    total_number_of_calls = 0
    for reader, writer in list(readers_to_writers.items()):
        variants_to_include = []

        headers = execution_context + reader.metaheaders
        for i, header in enumerate(headers):
            if not header.endswith("\n"):
                headers[i] += "\n"
        sorted_headers = utils.sort_metaheaders(headers)
        sorted_headers.append(reader.column_header + "\n")
        reader.open()

        for record in reader.tagged_vcf_records():
            if "JQ_EXCLUDE" not in record.filter:
                key = "^".join([record.chrom, record.pos, record.ref])
                if key in positions:
                    variants_to_include.append(record.text())

        _write_output(writer, sorted_headers, variants_to_include)
        total_number_of_calls += len(variants_to_include)

        reader.close()

    logger.info("Filtered to [{}] calls in high-confidence loci.",
                total_number_of_calls)
    logger.info("Jacquard wrote [{}] VCF files.",
                len(readers_to_writers.values()))

def _get_input_files(input_file):
    in_files = natsort.natsorted(glob.glob(os.path.join(input_file, "*.vcf")))
    if len(in_files) < 1:
        logger.error("Specified input directory [{}] contains no VCF files. "
                     "Check parameters and try again.", input_file)
        exit(1)

    logger.info("Processing [{}] VCF file(s) from [{}]",
                len(in_files),
                input_file)

    return in_files

def _write_output(writer, headers, variants):
    writer.open()
    for line in headers:
        writer.write(line)
    for line in variants:
        writer.write(line)
    writer.close()

def _build_readers(input_files):
    vcf_readers = []
    failures = 0
    for file_name in input_files:
        try:
            vcf_readers.append(vcf.VcfReader(vcf.FileReader(file_name)))
        except utils.JQException:
            failures += 1

    if failures:
        raise utils.JQException(("[{}] VCF files could not be parsed."
                                 " Review logs for details, adjust input, "
                                 "and try again.").format(failures))

    return vcf_readers

def _build_readers_to_writers(vcf_readers, output_dir):
    readers_to_writers = collections.OrderedDict()
    for reader in natsort.natsorted(vcf_readers):
        new_filename = _mangle_output_filenames(reader.file_name)
        output_filepath = os.path.join(output_dir, new_filename)

        readers_to_writers[reader] = vcf.FileWriter(output_filepath)

    return readers_to_writers

def _mangle_output_filenames(input_file):
    basename, extension = os.path.splitext(os.path.basename(input_file))
    return ".".join([basename, _FILE_OUTPUT_SUFFIX, extension.strip(".")])

def _get_output_filenames(input_files):
    output_files = set()
    for input_file in input_files:
        output_files.add(_mangle_output_filenames(input_file))

    return output_files

def _predict_output(args):
    input_file = os.path.abspath(args.input)
    input_files = sorted(glob.glob(os.path.join(input_file, "*.vcf")))

    desired_output_files = _get_output_filenames(input_files)

    return desired_output_files

def report_prediction(args):
    return _predict_output(args)

def get_required_input_output_types():
    return ("directory", "directory")

def add_subparser(subparser):
    # pylint: disable=line-too-long
    parser = subparser.add_parser("filter", help="Accepts a directory of Jacquard-tagged VCF results from one or more callers and creates a new directory of VCFs, where rows have been filtered to contain only positions that were called high-confidence somatic in any VCF.")
    parser.add_argument("input", help="Path to directory containing VCFs. All VCFs in this directory must have Jacquard-specific tags (see jacquard.py tag for more info")
    parser.add_argument("output", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")
    parser.add_argument("--include_variants", choices=["passed", "somatic"], help=("passed: Only include variants which passed their respective filter\n"
                                                                                   "somatic: Only include somatic variants"))

    parser.add_argument("--include_loci", choices=["any_passed", "all_passed", "any_soamtic", "all_somatic"], help=("any_passed: Include all variants at loci where at least one variant passed\n"
    # pylint: disable=line-too-long
                                                                                                                        "all_passed: Include all variants at loci where all variants passed\n"
                                                                                                                        "any_somatic: Include all variants at loci where at least one variant was somatic\n"
                                                                                                                        "all_somatic: Include all variants at loci where all variants were somatic"))
    parser.add_argument("--log_file", help="Log file destination")

def _validate_arguments(args):
    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)

    input_files = sorted(glob.glob(os.path.join(input_dir, "*.vcf")))
    out_files = sorted(glob.glob(os.path.join(output_dir, "*")))

    if len(input_files) < 1:
        logger.error("Specified input directory [{}] contains no VCF files. "
                     "Check parameters and try again.", input_dir)
        exit(1)

    full_path_input_files = [os.path.join(input_dir, i) for i in input_files]
    vcf_readers = _build_readers(full_path_input_files)
    writers_to_readers = _build_readers_to_writers(vcf_readers, output_dir)

    logger.info("Processing [{}] VCF file(s) from [{}]",
                len(input_files),
                input_dir)

    return writers_to_readers, out_files

def validate_args(dummy):
    pass

def execute(args, execution_context):
    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)

    include_variants = args.include_variants
    include_loci = args.include_loci

    in_files = _get_input_files(input_dir)
    file_readers = [vcf.FileReader(i) for i in in_files]
    factory = variant_caller_factory.VariantCallerFactory()
    dummy, trans_vcf_readers = factory.claim(file_readers)

    positions = _find_positions(trans_vcf_readers,
                                include_variants,
                                include_loci)

    execution_context = _add_headers(execution_context, positions)
    readers_to_writers = _build_readers_to_writers(trans_vcf_readers,
                                                   output_dir)

    _write_positions(readers_to_writers,
                     positions,
                     execution_context,
                     include_loci)

#pylint: disable=too-many-locals
from __future__ import print_function, absolute_import

from collections import defaultdict
import glob
import os
import re

import jacquard.logger as logger
import jacquard.utils as utils
import jacquard.variant_callers.variant_caller_factory as variant_caller_factory
import jacquard.vcf as vcf

JQ_OUTPUT_SUFFIX = "HCsomatic"
#TODO: (cgates): This module contains lots of file parsing/processing which should be using vcfReader structures
#TODO: (cgates): Use tuples instead of concatenated strings

#TODO: (cgates): refactor this as a stats object that collects info in the main processing loop
def _iterate_file(vcf_reader, num_records, somatic_positions, somatic):
    filtered_records = 0
    vcf_reader.open()
    for record in vcf_reader.vcf_records():
        if "JQ_EXCLUDE" in record.filter:
            filtered_records += 1
        else:
            sample_tag_values = record.sample_tag_values
            num_records += 1
            for sample in sample_tag_values:
                for tag in sample_tag_values[sample]:
                    if re.search(utils.jq_somatic_tag, tag):
                        if sample_tag_values[sample][tag] == "1":
                            somatic_key = "^".join([record.chrom,
                                                    record.pos,
                                                    record.ref])
                            somatic_positions[somatic_key] = 1
                            somatic = 1
    vcf_reader.close()

    return filtered_records, num_records, somatic_positions, somatic

def _find_somatic_positions(in_files):
    somatic_positions = {}
    no_jq_tags = []

    total_number_of_files = len(in_files)
    count = 1

    num_records = 0

    callers = defaultdict(int)
    for input_file in in_files:
        logger.info("Reading [{}] ({}/{})", os.path.basename(input_file), count,
                    total_number_of_files)
        somatic = 0
        vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
        caller = variant_caller_factory.get_caller(vcf_reader.metaheaders,
                                                   vcf_reader.column_header,
                                                   vcf_reader.file_name)

        (filtered_records,
         num_records,
         somatic_positions,
         somatic) =_iterate_file(vcf_reader,
                                 num_records,
                                 somatic_positions,
                                 somatic)

        callers[caller.name] += filtered_records

        if somatic == 0:
            no_jq_tags.append(input_file)
            logger.warning("Input file [{}] has no high-confidence somatic "
                           "variants.", os.path.basename(input_file))

#        in_file.close()

        count += 1

    total_filtered_records = 0

    for caller in callers:
        total_filtered_records += callers[caller]
        if callers[caller]:
            logger.debug("Removed [{}] problematic {} variant records with "
                         "filter=JQ_EXCLUDE", callers[caller], caller)
    if total_filtered_records:
        logger.warning("A total of [{}] problematic variant records failed "
                       "Jacquard's filters. See output and log for details.",
                       total_filtered_records)
    if no_jq_tags:
        logger.warning("[{}] VCF file(s) had no high-confidence somatic "
                       "variants. See log for details.", len(no_jq_tags))

    total_variants_searched = num_records + total_filtered_records
    logger.info("Searched [{}] variant calls.", total_variants_searched)

    logger.info("Found [{}] distinct high-confidence somatic positions.",
                len(somatic_positions.keys()))

    somatic_positions_header = "##jacquard.filterHCSomatic."\
                               "total_highConfidence_somatic_positions={0}\n"\
                               .format(len(somatic_positions.keys()))

    return somatic_positions, somatic_positions_header

def _write_somatic(in_files, output_file, somatic_positions, execution_context):
    total_number_of_calls = 0

    for input_file in in_files:
        non_somatic = 0
        headers = []
        actual_sorted_variants = []

        new_file = _mangle_output_filenames(input_file)

        in_file = open(input_file, "r")

        out_file = open(os.path.join(output_file, new_file), "w")

        for line in in_file:
            if line.startswith("#"):
                headers.append(line)
            elif "JQ_EXCLUDE" in line:
                continue
            else:
                split_line = line.split("\t")
                key = "^".join([split_line[0], split_line[1], split_line[3]])
                if key in somatic_positions:
                    actual_sorted_variants.append(line)
                else:
                    non_somatic += 1

        excluded_variants = "##jacquard.filterHCSomatic.excluded_variants="\
                            "{0}\n".format(non_somatic)
        headers.append(excluded_variants)
        headers.extend(execution_context)

        sorted_headers = _sort_headers(headers)
        for i, header in enumerate(sorted_headers):
            if not header.endswith("\n") and i != len(sorted_headers) - 1:
                sorted_headers[i] += "\n"

        _write_output(out_file, sorted_headers, actual_sorted_variants)
        total_number_of_calls += len(actual_sorted_variants)

        in_file.close()
        out_file.close()

    logger.info("Filtered to [{}] calls in high-confidence loci.",
                total_number_of_calls)

    logger.info("Jacquard wrote [{}] VCF files to [{}]",
                len(in_files),
                output_file)

def filter_somatic_positions(input_file, output_file, execution_context=None):
    if not execution_context:
        execution_context = []

    in_files = sorted(glob.glob(os.path.join(input_file, "*.vcf")))
    if len(in_files) < 1:
        logger.error("Specified input directory [{}] contains no VCF files. "
                     "Check parameters and try again.", input_file)
        exit(1)

    logger.info("Processing [{}] VCF file(s) from [{}]",
                len(in_files),
                input_file)

    (somatic_positions,
     somatic_positions_header) = _find_somatic_positions(in_files)
    execution_context.append(somatic_positions_header)

    _write_somatic(in_files, output_file, somatic_positions, execution_context)

def _write_output(writer, headers, actual_sorted_variants):
    for line in headers:
        writer.write(line)
    for line in actual_sorted_variants:
        writer.write(line)

def _sort_headers(headers):
    meta_headers = []
    field_header = ""
    for header in headers:
        if header.startswith("##"):
            header = header.replace("\t", "")
            meta_headers.append(header)
        else:
            field_header = header

    meta_headers.append(field_header)

    return meta_headers

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

def _build_writers_to_readers(vcf_readers, output_file):
    writers_to_readers = {}
    for reader in vcf_readers:
        new_filename = _mangle_output_filenames(reader.file_name)
        output_filepath = os.path.join(output_file, new_filename)

        writers_to_readers[reader] = vcf.FileWriter(output_filepath)

    return writers_to_readers

def _mangle_output_filenames(input_file):
    basename, extension = os.path.splitext(os.path.basename(input_file))
    return ".".join([basename, JQ_OUTPUT_SUFFIX, extension.strip(".")])
 
def _get_output_filenames(input_files):
    output_files = set()
    for input_file in input_files:
        output_files.add(_mangle_output_filenames(input_file))
 
    return output_files
 
def _predict_output(args):
    input_file = os.path.abspath(args.input)

    utils.validate_directories(input_file=input_file)
    input_files = sorted(glob.glob(os.path.join(input_file, "*.vcf")))

    desired_output_files = _get_output_filenames(input_files)

    return desired_output_files

def report_prediction(args):
    return _predict_output(args)

def add_subparser(subparser):
    # pylint: disable=line-too-long
    parser = subparser.add_parser("filter_hc_somatic", help="Accepts a directory of Jacquard-tagged VCF results from one or more callers and creates a new directory of VCFs, where rows have been filtered to contain only positions that were called high-confidence somatic in any VCF.")
    parser.add_argument("input", help="Path to directory containing VCFs. All VCFs in this directory must have Jacquard-specific tags (see jacquard.py tag for more info")
    parser.add_argument("output", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def _validate_arguments(args):
    input_file = os.path.abspath(args.input)
    output_file = os.path.abspath(args.output)

    utils.validate_directories(input_file, output_file)

    input_files = sorted(glob.glob(os.path.join(input_file, "*.vcf")))
    out_files = sorted(glob.glob(os.path.join(output_file, "*")))

    if len(input_files) < 1:
        logger.error("Specified input directory [{}] contains no VCF files. "
                     "Check parameters and try again.", input_file)
        exit(1)

    full_path_input_files = [os.path.join(input_file, i) for i in input_files]
    vcf_readers = _build_readers(full_path_input_files)
    writers_to_readers = _build_writers_to_readers(vcf_readers, output_file)

    logger.info("Processing [{}] VCF file(s) from [{}]",
                len(input_files),
                input_file)

    return writers_to_readers, out_files

def execute(args, execution_context):
    input_file = os.path.abspath(args.input)
    output_file = os.path.abspath(args.output)

    utils.validate_directories(input_file, output_file)

    filter_somatic_positions(input_file, output_file, execution_context)

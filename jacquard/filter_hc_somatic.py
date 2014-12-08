# pylint: disable=C0111
# pylint: disable-msg=W0403
import glob
import os
import re

import utils
import variant_callers.variant_caller_factory as variant_caller_factory
import vcf as vcf
import logger as logger
from collections import defaultdict

def _iterate_file(in_file, vcf_reader, filtered_records, num_records, somatic_positions, somatic):
    vcf_reader.open()
    for record in vcf_reader.vcf_records():
        if "JQ_EXCLUDE" in record.filter:
            filtered_records += 1
        else:
            num_records += 1
            for key in record.sample_dict:
                for tag in record.sample_dict[key]:
                    if re.search(utils.jq_somatic_tag, tag):
                        if record.sample_dict[key][tag] == "1":
                            somatic_key = "^".join([record.chrom,
                                                    record.pos])
                            somatic_positions[somatic_key] = 1
                            somatic = 1
    vcf_reader.close()

    return filtered_records, num_records, somatic_positions, somatic

def _find_somatic_positions(in_files, output_dir):
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

        in_file = open(input_file, "r")
        filtered_records = 0

        filtered_records, num_records, somatic_positions, somatic =_iterate_file(in_file, vcf_reader, filtered_records, num_records, somatic_positions, somatic)
        callers[caller.name] += filtered_records

        if somatic == 0:
            no_jq_tags.append(input_file)
            logger.warning("Input file [{}] has no high-confidence somatic "
                           "variants.", os.path.basename(input_file))

        in_file.close()

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

def _write_somatic(in_files, output_dir, somatic_positions, execution_context):
    total_number_of_calls = 0

    for input_file in in_files:
        non_somatic = 0
        headers = []
        actual_sorted_variants = []

        fname, extension = os.path.splitext(os.path.basename(input_file))
        new_file = fname + "_HCsomatic" + extension
        in_file = open(input_file, "r")

        out_file = open(os.path.join(output_dir, new_file), "w")

        for line in in_file:
            if line.startswith("#"):
                headers.append(line)
            elif "JQ_EXCLUDE" in line:
                continue
            else:
                split_line = line.split("\t")
                key = "^".join([split_line[0], split_line[1]])
                if key in somatic_positions:
                    actual_sorted_variants.append(line)
                else:
                    non_somatic += 1

        excluded_variants = "##jacquard.filterHCSomatic.excluded_variants="\
                            "{0}\n".format(non_somatic)
        headers.append(excluded_variants)
        headers.extend(execution_context)

        sorted_headers = utils.sort_headers(headers)
        for i,header in enumerate(sorted_headers):
            if not header.endswith("\n") and i!=len(sorted_headers)-1:
                sorted_headers[i]+="\n"

        utils.write_output(out_file, sorted_headers, actual_sorted_variants)
        total_number_of_calls += len(actual_sorted_variants)

        in_file.close()
        out_file.close()

    logger.info("Filtered to [{}] calls in high-confidence loci.",
                total_number_of_calls)

    logger.info("Jacquard wrote [{}] VCF files to [{}]",
                len(in_files),
                output_dir)

def filter_somatic_positions(input_dir, output_dir, execution_context=[]):
    in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
    if len(in_files) < 1:
        logger.error("Specified input directory [{}] contains no VCF files. "
                     "Check parameters and try again.", input_dir)
        exit(1)

    logger.info("Processing [{}] VCF file(s) from [{}]",
                len(in_files),
                input_dir)

# pylint: disable=C0301
    somatic_positions, somatic_positions_header = _find_somatic_positions(in_files,
                                                                          output_dir)
    execution_context.append(somatic_positions_header)

    _write_somatic(in_files, output_dir, somatic_positions, execution_context)


def add_subparser(subparser):
    # pylint: disable=C0301
    parser = subparser.add_parser("filter_hc_somatic", help="Accepts a directory of Jacquard-tagged VCF results from one or more callers and creates a new directory of VCFs, where rows have been filtered to contain only positions that were called high-confidence somatic in any VCF.")
    parser.add_argument("input", help="Path to directory containing VCFs. All VCFs in this directory must have Jacquard-specific tags (see jacquard.py tag for more info")
    parser.add_argument("output", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)

    utils.validate_directories(input_dir, output_dir)

    filter_somatic_positions(input_dir, output_dir, execution_context)

# pylint: disable=W0212,C0111
from collections import OrderedDict
import glob
import os
import re

import logger as logger
import utils as utils
import vcf as vcf

UNUSED_REGEX_WARNING_FORMAT = ("The expression [{}] in column specification " +
                               "file [{}:{}] didn't match any input columns; " +
                               "columns may have matched earlier expressions, "+
                               "or this expression may be irrelevant.")

def _read_col_spec(col_spec):
    if not os.path.isfile(col_spec):
        raise utils.JQException("The column specification file [{}] could "
                                "not be read. "
                                "Review inputs/usage and try again.",
                                col_spec)

    spec_file = open(col_spec, "r")
    columns = []

    for line in spec_file:
        columns.append(line.strip("\n"))

    spec_file.close()

    return columns

def _path_type(path):
    return "file" if os.path.isfile(path) else "directory"

def _build_output_file_names(input_path, output_path):
    input_files = sorted(glob.glob(os.path.join(input_path, "*.vcf")))
    if len(input_files) == 0:
        raise utils.JQException(("Specified input directory [{}] contains "
                                 "no VCF files. Review inputs and try "
                                 "again."),
                                input_path)

    basenames = [os.path.splitext(os.path.basename(i))[0] + ".txt" \
                       for i in input_files]
    output_path = [os.path.join(output_path, i) for i in basenames]

    return input_files, output_path


def _validate_input_and_output(input_path, output_path):
    input_path_type = _path_type(input_path)
    if os.path.exists(output_path) and \
            input_path_type != _path_type(output_path):
        raise utils.JQException(("Specified output [{0}] must be a {1} "
                                 "if input [{2}] is a {1}. Review "
                                 "arguments and try again."),
                                output_path,
                                input_path_type,
                                input_path)
    if os.path.isfile(input_path):
        return [input_path], [output_path]
    else:
        return _build_output_file_names(input_path, output_path)

##TODO: hook this idea up -- change method
def _disambiguate_column_names(column_header, info_header):
    overlap = 0
    for column in info_header:
        if column in column_header:
            overlap = 1

    return ["INFO_" + i for i in info_header] if overlap else info_header

def _create_row_dict(column_list, vcf_record):
    row_dict = {"CHROM" : vcf_record.chrom,
                "POS" : vcf_record.pos,
                "ID" : vcf_record.id,
                "REF" : vcf_record.ref,
                "ALT" : vcf_record.alt,
                "QUAL" : vcf_record.qual,
                "FILTER" : vcf_record.filter}

    for i, sample_name in enumerate(column_list[9:]):
        format_key_values = vcf_record.sample_dict[i]
        for format_key, format_value in format_key_values.items():
            row_dict[format_key + "|" + sample_name] = format_value

    row_dict = dict(row_dict.items() + vcf_record.get_info_dict().items())

    return row_dict

def _create_actual_column_list(column_spec_list,
                               potential_col_list,
                               column_spec_filename):

    actual_column_list = []
    for i, column_regex in enumerate(column_spec_list):
        no_columns_found = True
        for column_name in potential_col_list:
            if re.match(column_regex, column_name):
                actual_column_list.append(column_name)
                no_columns_found = False

        if no_columns_found:
            logger.warning(UNUSED_REGEX_WARNING_FORMAT,
                           column_regex,
                           column_spec_filename,
                           i + 1
                           )

    if actual_column_list:
        return actual_column_list
    else:
        raise utils.JQException("The column specification file [{}] would "
                                "exclude all input columns. Review "
                                "inputs/usage and try again.",
                                column_spec_filename)

def _create_potential_column_list(vcf_reader):
    column_headers = vcf_reader.get_col_header_list()
    info_dict = vcf_reader.info_metaheaders
    format_dict = vcf_reader.format_metaheaders

    format_sample_names = []
    for sample_name in column_headers[9:]:
        for format_tag in format_dict.keys():
            format_sample_names.append(format_tag + "|" + sample_name)

    static_column_headers = column_headers[0:10]
    processed_column_headers = static_column_headers + format_sample_names

    return processed_column_headers + info_dict.keys() + format_dict.keys()

# pylint: disable=C0301
def add_subparser(subparser):
    parser = subparser.add_parser("expand", help="Pivots annotated VCF file so that given sample specific information is fielded out into separate columns. Returns an Excel file containing concatenation of all input files.")
    parser.add_argument("input", help="Path to annotated VCF file or path to directory of annotated VCF files. Other file types ignored")
    parser.add_argument("output", help="Path to directory of output variant-level TXT files")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("-c", "--column_specification", help="Path to text file containing column regular expressions to be included in output file")
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def execute(args, execution_context):
    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)
    col_spec = args.column_specification if args.column_specification else 0

    col_spec_columns = _read_col_spec(col_spec) if col_spec else 0
    input_files, output_files = _validate_input_and_output(input_path,
                                                           output_path)

    logger.info("Expanding {} VCF files in [{}] to [{}]",
                len(input_files),
                input_path,
                output_path)

    for i, input_file in enumerate(input_files):
        output_file = output_files[i]
        file_reader = vcf.FileReader(input_file)
        vcf_reader = vcf.VcfReader(file_reader)

        file_writer = vcf.FileWriter(output_file)
        file_writer.open()

        potential_columns = _create_potential_column_list(vcf_reader)

        if col_spec_columns:
            actual_columns = _create_actual_column_list(col_spec_columns,
                                                        potential_columns,
                                                        col_spec)
        else:
            actual_columns = potential_columns

        file_writer.write("#" + "\t".join(actual_columns) + "\n")

        vcf_reader.open()
        for vcf_record in vcf_reader.vcf_records():
            original_col_header = vcf_reader.get_col_header_list()
            row_dict = _create_row_dict(original_col_header, vcf_record)

            new_line = []
            for col in actual_columns:
                if col in row_dict:
                    new_line.append(row_dict[col])
                else:
                    new_line.append(".")

            file_writer.write("\t".join(new_line) + "\n")
        file_writer.close()
        vcf_reader.close()

    logger.info("Wrote [{}] VCF files to [{}]", len(input_files), output_path)


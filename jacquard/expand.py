#pylint: disable=unused-argument, too-many-locals
from __future__ import print_function, absolute_import
import os
import re

import jacquard.logger as logger
import jacquard.utils as utils
import jacquard.vcf as vcf

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
        columns.append(line.rstrip())

    spec_file.close()

    return columns


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

    for sample_name in column_list[9:]:
        format_key_values = vcf_record.sample_tag_values[sample_name]
        for format_key, format_value in format_key_values.items():
            row_dict[format_key + "|" + sample_name] = format_value

    row_dict = dict(row_dict.items() + vcf_record.info_dict.items())

    return row_dict

def _create_actual_column_list(column_spec_list,
                               potential_col_list,
                               column_spec_filename):

    actual_column_list = []
    for i, column_regex in enumerate(column_spec_list):
        no_columns_found = True
        for column_name in potential_col_list:
            column_exists = column_name in actual_column_list
            regex_matches = re.match("^" + column_regex + "$", column_name)
            if regex_matches and not column_exists:
                actual_column_list.append(column_name)
                no_columns_found = False

        if no_columns_found:
            logger.warning(UNUSED_REGEX_WARNING_FORMAT,
                           column_regex,
                           column_spec_filename,
                           i + 1)

    if actual_column_list:
        return actual_column_list
    else:
        raise utils.JQException("The column specification file [{}] would "
                                "exclude all input columns. Review "
                                "inputs/usage and try again.",
                                column_spec_filename)

def _create_potential_column_list(vcf_reader):
    format_sample_names = []
    for format_tag in sorted(vcf_reader.format_metaheaders.keys()):
        for sample_name in vcf_reader.sample_names:
            format_sample_names.append(format_tag + "|" + sample_name)

    static_column_headers = vcf_reader.split_column_header[0:8]

    return static_column_headers \
           + sorted(vcf_reader.info_metaheaders.keys()) \
           + format_sample_names

def add_subparser(subparser):
    # pylint: disable=C0301
    parser = subparser.add_parser("expand", help="Pivots annotated VCF file so that given sample specific information is fielded out into separate columns. Returns an Excel file containing concatenation of all input files.")
    parser.add_argument("input", help="Path to annotated VCF file or path to directory of annotated VCF files. Other file types ignored")
    parser.add_argument("output", help="Path to directory of output variant-level TXT files")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("-c", "--column_specification", help="Path to text file containing column regular expressions to be included in output file")
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def _predict_output(args):
    return set([os.path.basename(args.output)])

def report_prediction(args):
    return _predict_output(args)

def get_required_input_output_types():
    return ("file", "file")

def execute(args, execution_context):
    input_file = os.path.abspath(args.input)
    output_file = os.path.abspath(args.output)

    col_spec = args.column_specification if args.column_specification else 0

    col_spec_columns = _read_col_spec(col_spec) if col_spec else 0

    logger.info("Expanding [{}] to [{}]",
                input_file,
                output_file)

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
        row_dict = _create_row_dict(vcf_reader.split_column_header,
                                    vcf_record)

        new_line = []
        for col in actual_columns:
            if col in row_dict:
                new_line.append(row_dict[col])
            else:
                new_line.append(".")

        file_writer.write("\t".join(new_line) + "\n")
    file_writer.close()
    vcf_reader.close()

    logger.info("Wrote input [{}] to output [{}]", input_file, output_file)

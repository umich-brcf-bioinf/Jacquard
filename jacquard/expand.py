"""Expands a VCF file to tab-separated text file

The number of rows will be the same as the number of records in the VCF (i.e.
line count of VCF including the column header - number of metaheaders).

Will create a column for each:
 * static field (e.g. CHROM, POS, ID, REF, ALT, QUAL, FILTER)
 * INFO field listed in the metaheader
 * combination of sample name and FORMAT tag listed in the metaheader

Expand relies on accurate metaheaders; it will not expand any fields absent
from the metaheaders.
"""
from __future__ import print_function, absolute_import, division

import argparse
from collections import OrderedDict
import os
import re

import jacquard.utils.logger as logger
import jacquard.utils.utils as utils
import jacquard.utils.vcf as vcf


UNUSED_REGEX_WARNING_FORMAT = ("The expression [{}] in column specification "
                               "file [{}:{}] didn't match any input columns; "
                               "columns may have matched earlier expressions, "
                               "or this expression may be irrelevant.")

def _read_col_spec(col_spec):
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

#TODO: cgates: suspect this could be simplified?
def _create_row_dict(column_list, vcf_record):
    row_dict = {"CHROM" : vcf_record.chrom,
                "POS" : vcf_record.pos,
                "ID" : vcf_record.vcf_id,
                "REF" : vcf_record.ref,
                "ALT" : vcf_record.alt,
                "QUAL" : vcf_record.qual,
                "FILTER" : vcf_record.filter}

    for sample_name in column_list[9:]:
        format_key_values = vcf_record.sample_tag_values[sample_name]
        for format_key, format_value in format_key_values.items():
            row_dict[format_key + "|" + sample_name] = format_value

    new_dict = row_dict.copy()
    new_dict.update(vcf_record.info_dict)

    return new_dict

def _filter_column_list(column_spec_list,
                        potential_col_list,
                        column_spec_filename):

    actual_column_list = []
    glossary_fields = OrderedDict()
    for i, column_regex in enumerate(column_spec_list):
        no_columns_found = True
        for (column_name, vcf_tag) in potential_col_list.items():
            column_exists = column_name in actual_column_list
            regex_matches = re.match("^" + column_regex + "$", column_name)
            if regex_matches and not column_exists:
                glossary_fields[vcf_tag] = 1
                actual_column_list.append(column_name)
                no_columns_found = False

        if no_columns_found:
            logger.warning(UNUSED_REGEX_WARNING_FORMAT,
                           column_regex,
                           column_spec_filename,
                           i + 1)

    if not actual_column_list:
        raise utils.JQException("The column specification file [{}] would "
                                "exclude all input columns. Review "
                                "inputs/usage and try again.",
                                column_spec_filename)
    glossary_fields = [x for x in glossary_fields.keys() if x]
    return actual_column_list, glossary_fields


def _create_potential_column_list(vcf_reader):
    fixed_fields = vcf_reader.split_column_header[0:8]
    columns = OrderedDict(zip(fixed_fields, [None] * 8))

    for info_tag in sorted(vcf_reader.info_metaheaders.keys()):
        columns[info_tag] = info_tag

    for format_tag in sorted(vcf_reader.format_metaheaders.keys()):
        for sample_name in vcf_reader.sample_names:
            columns[format_tag + "|" + sample_name] = format_tag

    return columns

def _get_glossary_writer(output_file):
    output_dir = os.path.dirname(output_file)
    complete_output_fname = os.path.basename(output_file)
    output_fname, dummy = os.path.splitext(complete_output_fname)
    glossary_fname = ".".join([output_fname, "glossary.txt"])
    glossary = os.path.join(output_dir, glossary_fname)

    return vcf.FileWriter(glossary)

def _create_glossary(metaheaders, glossary_fields, writer):
    glossary_entries = dict([_create_glossary_entry(x) for x in metaheaders])
    writer.write("FIELD_NAME\tTYPE\tDESCRIPTION\n")
    for field_name in glossary_fields:
        writer.write(glossary_entries[field_name])

def _create_glossary_entry(metaheader):
    try:
        header_type = re.search(r'^##(.*)=<', metaheader).group(1)
        id_value = re.search(r'^##.*\WID=(\w*)', metaheader).group(1)
        desc = re.search(r'^##.*Description="([^"]*)', metaheader).group(1)
        entry = "\t".join([id_value, header_type, desc.strip('"')]) + "\n"
        return (id_value, entry)
    except AttributeError:
        return (None, None)


def add_subparser(subparser):
    # pylint: disable=C0301
    parser = subparser.add_parser("expand", formatter_class=argparse.RawTextHelpFormatter, help="Pivots annotated VCF file so that given sample specific information is fielded out into separate columns. Returns an Excel file containing concatenation of all input files.")
    parser.add_argument("input", help="Path to annotated VCF file or path to directory of annotated VCF files. Other file types ignored")
    parser.add_argument("output", help="Path to directory of output variant-level TXT files")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("-c", "--column_specification",
                        help="Path to text file containing column regular expressions to be included in output file",)
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")
    parser.add_argument("--log_file", help="Log file destination")

def _predict_output(args):
    return set([os.path.basename(args.output)])

def report_prediction(args):
    return _predict_output(args)

def get_required_input_output_types():
    return ("file", "file")

def validate_args(args):
    if args.column_specification:
        if not os.path.isfile(args.column_specification):
            raise utils.UsageError(("The column specification file [{}] could "
                                    "not be read. Review inputs/usage and "
                                    "try again."),
                                    args.column_specification)

        _read_col_spec(args.column_specification)

def _get_actual_columns(vcf_reader, col_spec):
    columns  = _create_potential_column_list(vcf_reader)
    glossary_fields = []
    for glossary_field in sorted([x for x in columns.values() if x], 
                                 key=lambda x: x.upper()):
        if glossary_field and glossary_field not in glossary_fields:
            glossary_fields.append(glossary_field)

    if col_spec:
        col_spec_columns = _read_col_spec(col_spec)
        (columns,
         glossary_fields) = _filter_column_list(col_spec_columns,
                                                columns,
                                                col_spec)
    return columns, glossary_fields


def execute(args, dummy_execution_context):
    #for the moment, there is no good place to put the execution context
    input_file = os.path.abspath(args.input)
    output_file = os.path.abspath(args.output)
    col_spec = args.column_specification if args.column_specification else None

    logger.debug("Expanding [{}] to [{}]",
                 input_file,
                 output_file)
    logger.info("Expanding [{}] to [{}]",
                args.input,
                args.original_output)

    vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
    file_writer = vcf.FileWriter(output_file)
    file_writer.open()

    (columns, glossary_fields) = _get_actual_columns(vcf_reader, col_spec)

    file_writer.write("#" + "\t".join(columns) + "\n")

    line_count = 0
    vcf_reader.open()
    for vcf_record in vcf_reader.vcf_records():
        row_dict = _create_row_dict(vcf_reader.split_column_header,
                                    vcf_record)

        new_line = []
        for col in columns:
            if col in row_dict:
                new_line.append(row_dict[col])
            else:
                new_line.append(".")

        file_writer.write("\t".join(new_line) + "\n")
        line_count +=1
        if line_count % 10000 == 0:
            logger.info("Expanding: {} rows processed", line_count)
    logger.info("Expand complete: {} rows processed", line_count)

    file_writer.close()

    glossary_writer = _get_glossary_writer(output_file)
    glossary_writer.open()
    _create_glossary(vcf_reader.metaheaders, glossary_fields, glossary_writer)
    glossary_writer.close()
    logger.info("Wrote glossary to [{}]",
                os.path.basename(glossary_writer.output_filepath))

    vcf_reader.close()
    logger.debug("Wrote input [{}] to output [{}]", input_file, output_file)

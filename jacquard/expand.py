from collections import OrderedDict
import glob
import os
import re

import logger as logger
import utils as utils
import vcf as vcf

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

def _validate_input_and_output(input_path, output_path):
    ##input is a file
    if os.path.isfile(input_path):
        if os.path.isdir(output_path):
            raise utils.JQException("Specified output {} must be a file "
                                    "if input {} is given as a file. Review "
                                    "inputs and try again.", output_path,
                                    input_path)
        return [input_path], [output_path]

##input is a directory
    elif os.path.isdir(input_path):
        input_files= sorted(glob.glob(os.path.join(input_path,"*.vcf")))
        if len(input_files) == 0:
            raise utils.JQException("Specified input directory {} contains "+
                                    "no VCF files. Review inputs and try again.",
                                    input_path)

        if os.path.isfile(output_path):
            raise utils.JQException("Specified output {} must be a directory "
                                    "if input {} is given as a directory."
                                    "Review inputs and try again.", output_path,
                                    input_path)
        try:
            os.mkdir(output_path)
        except:
            pass

        tmp_output_path = [os.path.splitext(os.path.basename(i))[0] + ".txt" for i in input_files]
        output_path = [os.path.join(output_path, i) for i in tmp_output_path]

        return input_files, output_path

def _parse_meta_headers(meta_headers):
    info_fields = []
    format_tags = []

    for meta_header in meta_headers:
        if meta_header.startswith("##FORMAT="):
            tag = meta_header.split(",")[0].split("=")[-1]
            format_tags.append(tag)
        elif meta_header.startswith("##INFO="):
            field = meta_header.split(",")[0].split("=")[-1]
            info_fields.append(field)

    info_fields.sort()
    format_tags.sort()

    if len(info_fields) == 0 or len(format_tags) == 0:
        raise utils.JQException("Unable to parse meta_headers for INFO and/" +
                                "or FORMAT fields. Review input and try again.")

    return info_fields, format_tags

def _append_format_tags_to_samples(format_tags, samples):
    format_sample_header = []

    for sample in samples:
        for tag in format_tags:
            modified_field = tag + "|" + sample
            format_sample_header.append(modified_field)

    return format_sample_header

def _create_filtered_header(header_dict, filtered_header_dict, desired_column):
    not_found = 1

    for key in header_dict.keys():
        for incoming_column in header_dict[key]:
            if re.search(desired_column, incoming_column):
                if key in filtered_header_dict:
                    filtered_header_dict[key].append(incoming_column)
                else:
                    filtered_header_dict[key] = [incoming_column]

                not_found = 0

    return filtered_header_dict, not_found

def _validate_column_specification(filtered_header_dict, not_found_regex):
    invalid = 1
    for val in filtered_header_dict.values():
        if len(val) != 0:
            invalid = 0

    if invalid:
        raise utils.JQException("The column specification file would "
                                "exclude all input columns. Review inputs/"
                                "usage and try again.")

    if len(not_found_regex) != 0:
        logger.warning("The expression {} in column specification file didn't "
                       "match any input columns. Columns may have matched "
                       "earlier expressions, or this expression may be "
                       "irrelevant.", not_found_regex)

def _filter_and_sort(header_dict, columns_to_expand):
    filtered_header_dict = OrderedDict()
    not_found_regex = []

    for desired_column in columns_to_expand:
        filtered_header_dict, not_found = _create_filtered_header(header_dict,
                                                       filtered_header_dict,
                                                       desired_column)
        if not_found:
            not_found_regex.append(desired_column)

    _validate_column_specification(filtered_header_dict, not_found_regex)

    return filtered_header_dict

def _get_headers(vcf_reader):
    split_column_header = vcf_reader.column_header.split("\t")

    column_header_no_samples = split_column_header[0:7]
    samples = split_column_header[9:]

    (info_header, format_tags) = _parse_meta_headers(vcf_reader.metaheaders)
    format_header = _append_format_tags_to_samples(format_tags, samples)

    return column_header_no_samples, info_header, format_header

def _disambiguate_column_names(column_header, info_header):
    overlap = 0
    for column in info_header:
        if column in column_header:
            overlap = 1

    return ["INFO_" + i for i in info_header] if overlap else info_header

def _parse_info_field(vcf_record, info_header):
    info_dict = vcf_record.get_info_dict()
    info_columns = []

    for tag in info_header:
        info_cell = ""
        if tag in info_dict:
            info_cell = info_dict[tag]

        info_columns.append(info_cell)

    return info_columns

def _parse_format_tags(vcf_record, format_header, column_header):
    samples = column_header.split("\t")[9:]
    sample_dict = vcf_record.sample_dict
    format_columns = []

    for key in sample_dict.keys():
        sample_dict[samples[key]] = sample_dict.pop(key)

    for field in format_header:
        tag = field.split("|")[0]
        sample_name = "|".join(field.split("|")[1:])

        try:
            format_columns.append(sample_dict[sample_name][tag])
        except KeyError:
            format_columns.append("")

    return format_columns

def _write_vcf_records(vcf_reader, file_writer, header_dict):
    vcf_reader.open()

    for record in vcf_reader.vcf_records():
        row = []
        for key in header_dict.keys():
            if key == "column_header":
                column_header = [i.strip("#").lower() for i in header_dict[key]]
                column_values = [getattr(record, i) for i in column_header]
                row.extend(column_values)

            elif key == "info_header":
                info_columns = _parse_info_field(record, header_dict[key])
                row.extend(info_columns)

            elif key == "format_header":
                format_columns = _parse_format_tags(record, header_dict[key],
                                            vcf_reader.column_header)
                row.extend(format_columns)

        row_string = "\t".join(row) + "\n"
        file_writer.write(row_string)

    vcf_reader.close()

def _create_complete_header(header_dict):
    complete_header = []
    for i in header_dict.values():
        complete_header.extend(i)

    return "\t".join(complete_header) + "\n"

# pylint: disable=C0301
def add_subparser(subparser):
    parser = subparser.add_parser("expand", help="Pivots annotated VCF file so that given sample specific information is fielded out into separate columns. Returns an Excel file containing concatenation of all input files.")
    parser.add_argument("input", help="Path to annotated VCF file or path to directory of annotated VCF files. Other file types ignored")
    parser.add_argument("output", help="Path to directory of output variant-level TXT files")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("-c", "--column_specification", help="Path to text file containing column regular expressions to be included in output file")

#     return parser

def execute(args, execution_context):
    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)
    col_spec = args.column_specification if args.column_specification else 0

    columns_to_expand = _read_col_spec(col_spec) if col_spec else 0

    input_files, output_files = _validate_input_and_output(input_path, output_path)
    logger.info("Expanding {} VCF files in {} to {}", len(input_files), input_path, output_path)

    for i, input_file in enumerate(input_files):
        output_file = output_files[i]
        file_reader = vcf.FileReader(input_file)
        vcf_reader = vcf.VcfReader(file_reader)

        column_header, info_header, format_header = _get_headers(vcf_reader)

        info_header = _disambiguate_column_names(column_header, info_header)
        header_dict = OrderedDict([("column_header", column_header), 
                                   ("info_header", info_header),
                                   ("format_header", format_header)])
 
        if columns_to_expand:
            header_dict = _filter_and_sort(header_dict, columns_to_expand)

        file_writer = vcf.FileWriter(output_file)
        file_writer.open()

        logger.info("Writing {} to {}", input_file, output_file)
        file_writer.write(_create_complete_header(header_dict))
        _write_vcf_records(vcf_reader, file_writer, header_dict)

        file_writer.close()

    logger.info("Wrote {} VCF files to {}", len(input_files), output_path)


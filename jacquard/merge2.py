# pylint: disable=C0111
from collections import defaultdict, OrderedDict
import glob
import os
import re

import utils as utils
import vcf as vcf

def _produce_merged_metaheaders(vcf_reader, all_meta_headers, count):
    vcf_reader.open()
    for meta_header in vcf_reader.metaheaders:
        if meta_header not in all_meta_headers:
            all_meta_headers.append(meta_header)

    samples = vcf_reader.column_header.split("\t")[9:]
    all_meta_headers.append("##jacquard.merge.file{0}={1}({2})"\
                            .format(count, vcf_reader.file_name, samples))
    all_meta_headers.append('##INFO=<ID=JQ_MULT_ALT_LOCUS,Number=0,Type=Flag,'\
                            'Description="dbSNP Membership",Source="Jacquard",'\
                            'Version="{}">'.format(utils.__version__))
    vcf_reader.close()

    return all_meta_headers

def _get_coordinate_dict(vcf_reader, coordinate_dict):
    vcf_reader.open()

    for vcf_record in vcf_reader.vcf_records():
        key = "^".join([vcf_record.chrom,
                        vcf_record.pos])
        coordinate = "^".join([vcf_record.chrom,
                               vcf_record.pos,
                               vcf_record.ref,
                               vcf_record.alt])
        coordinate_dict[key].append(coordinate)

    vcf_reader.close()

    return coordinate_dict

def _sort_coordinate_dict(coordinate_dict):
    def _convert(element):
        return int(element) if element.isdigit() else element

    def _alphanum_key (coordinates):
        return [_convert(c) for c in re.split('([0-9]+)', coordinates[0])]

    return OrderedDict(sorted(coordinate_dict.items(), key=_alphanum_key))

def _write_metaheaders(file_writer,
                      all_metaheaders,
                      column_header,
                      execution_context):

    all_metaheaders.extend(execution_context)
    file_writer.write("\n".join(all_metaheaders) + "\n")
    file_writer.write("\t".join(column_header) + "\n")

def _alter_record_fields(vcf_record, coordinates):
    vcf_record.id = "."
    vcf_record.qual = "."
    vcf_record.filter = "."

    if len(coordinates) > 1:
        info = vcf_record.info.split(";") if vcf_record.info != "." else []
        info.append("JQ_MULT_ALT_LOCUS")
        vcf_record.info = ";".join(info)

    return vcf_record

def _write_variants(vcf_reader, file_writer, coordinates):
    vcf_reader.open()

    for vcf_record in vcf_reader.vcf_records():
        vcf_record = _alter_record_fields(vcf_record, coordinates)

        record_key = "^".join([vcf_record.chrom,
                        vcf_record.pos,
                        vcf_record.ref,
                        vcf_record.alt])

        for coordinate in coordinates:
            if record_key == coordinate:
                ##when we add FORMAT and SAMPLE fields, we don't have to call
                ##vcf_record.asText() with any arguments.
                file_writer.write(vcf_record.asText(stringifier=[vcf_record.chrom,
                                                                 vcf_record.pos,
                                                                 vcf_record.id,
                                                                 vcf_record.ref,
                                                                 vcf_record.alt,
                                                                 vcf_record.qual,
                                                                 vcf_record.filter,
                                                                 vcf_record.info]))
        else:
            continue

    vcf_reader.close()

def add_subparser(subparser):
    #pylint: disable=C0301
    parser = subparser.add_parser("merge", help="Accepts a directory of VCFs and returns a single merged VCF file.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to output variant-level VCF file")
    parser.add_argument("-a", "--allow_inconsistent_sample_sets", action="store_true", default=False, help="Allow inconsistent sample sets across callers. Not recommended.")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def execute(args, execution_context):
    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)

    all_metaheaders = []
    coordinate_dict = defaultdict(list)
    input_files = sorted(glob.glob(os.path.join(input_path, "*.vcf")))

    file_writer = vcf.FileWriter(output_path)
    file_writer.open()

    count = 1
    for input_file in input_files:
        vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
        all_metaheaders = _produce_merged_metaheaders(vcf_reader,
                                                      all_metaheaders,
                                                      count)

        column_header = vcf_reader.column_header.split("\t")[0:9]
        coordinate_dict = _get_coordinate_dict(vcf_reader, coordinate_dict)
        count += 1

    _write_metaheaders(file_writer,
                      all_metaheaders,
                      column_header,
                      execution_context)

    sorted_coordinates = _sort_coordinate_dict(coordinate_dict)
    for coordinates in sorted_coordinates.values():
        for input_file in input_files:
            vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
            _write_variants(vcf_reader, file_writer, coordinates)

    file_writer.close()



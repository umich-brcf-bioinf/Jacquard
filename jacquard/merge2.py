# pylint: disable=C0111
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

def _get_all_coordinates(vcf_reader, all_coordinates):
    vcf_reader.open()

    for vcf_record in vcf_reader.vcf_records():
        coordinate = "^".join([vcf_record.chrom,
                               vcf_record.pos,
                               vcf_record.ref,
                               vcf_record.alt])
        all_coordinates.append(coordinate)

    vcf_reader.close()

    return all_coordinates

def _sort_all_coordinates(all_coordinates):
    def _convert(element):
        return int(element) if element.isdigit() else element

    def _alphanum_key (coordinates):
        return [_convert(c) for c in re.split('([0-9]+)', coordinates)]

    return sorted(all_coordinates, key=_alphanum_key)

def _write_metaheaders(file_writer,
                      all_metaheaders,
                      column_header,
                      execution_context):

    all_metaheaders.extend(execution_context)
    file_writer.write("\n".join(all_metaheaders) + "\n")
    file_writer.write("\t".join(column_header) + "\n")

def _write_variants(vcf_reader, file_writer, coordinate):
    vcf_reader.open()

    for vcf_record in vcf_reader.vcf_records():
        vcf_record.id = "."
        vcf_record.qual = "."
        vcf_record.filter = "."

        key = "^".join([vcf_record.chrom,
                        vcf_record.pos,
                        vcf_record.ref,
                        vcf_record.alt])

        if key == coordinate:
            file_writer.write(vcf_record.asText())
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
    all_coordinates = []
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
        all_coordinates = _get_all_coordinates(vcf_reader, all_coordinates)
        count += 1

    _write_metaheaders(file_writer,
                      all_metaheaders,
                      column_header,
                      execution_context)

    sorted_coordinates = _sort_all_coordinates(all_coordinates)
    for coordinate in sorted_coordinates:
        for input_file in input_files:
            vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
            _write_variants(vcf_reader, file_writer, coordinate)

    file_writer.close()



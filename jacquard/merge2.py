# pylint: disable=C0111
from collections import defaultdict, OrderedDict
import glob
import os
import re
from sets import Set

import utils as utils
import vcf as vcf

class BufferedReader(object):
    def __init__(self, reader):
        self.reader = reader
        self.base_iterator = iter(reader.vcf_records())
        self.current_record = self.base_iterator.next()
        self.sample_dict = {}

    def check_current_value(self, val):
        if val == self.current_record:
            sample_names = self.reader.column_header.split("\t")[9:]
            for i, sample_name in enumerate(sample_names):
                key = self.reader.file_name + "|" + sample_name
                self.sample_dict[key] = self.current_record.sample_dict[i]
            self.current_record = self._get_next()

        return self.sample_dict

    def _get_next(self):
        try:
            return self.base_iterator.next()
        except StopIteration:
            return None

def _produce_merged_metaheaders(vcf_reader, all_meta_headers, count):
    vcf_reader.open()
    for meta_header in vcf_reader.metaheaders:
        if meta_header not in all_meta_headers:
#             if re.search('##FORMAT=.*Source="Jacquard")', meta_header):
            all_meta_headers.append(meta_header)

    samples = vcf_reader.column_header.split("\t")[9:]
    all_meta_headers.append("##jacquard.merge.file{0}={1}({2})"\
                            .format(count, vcf_reader.file_name, samples))
    all_meta_headers.append('##INFO=<ID=JQ_MULT_ALT_LOCUS,Number=0,Type=Flag,'\
                            'Description="dbSNP Membership",Source="Jacquard",'\
                            'Version="{}">'.format(utils.__version__))
    vcf_reader.close()

    return all_meta_headers

def _extract_format_ids(all_meta_headers):
    pattern = re.compile("##FORMAT=.*[<,]ID=([A-Za-z0-9_]+)")
    format_tags = Set()
    for meta_header in all_meta_headers:
        match = pattern.search(meta_header)
        if match:
            format_tags.add(match.group(1))

    return format_tags

def _add_to_coordinate_set(vcf_reader, coordinate_set):
    vcf_reader.open()
    for vcf_record in vcf_reader.vcf_records():
        coordinate_set.add(vcf_record.get_empty_record())

    vcf_reader.close()

    return coordinate_set

def _sort_coordinate_set(coordinate_set):
    coordinate_list = list(coordinate_set)
    coordinate_list.sort()
    return coordinate_list
     
def _write_metaheaders(file_writer,
                      all_metaheaders,
                      column_header,
                      execution_context):

    all_metaheaders.extend(execution_context)
    file_writer.write("\n".join(all_metaheaders) + "\n")
    file_writer.write("\t".join(column_header) + "\n")

def _create_reader_lists(input_files):
    buffered_readers = []
    vcf_readers = []

    for input_file in input_files:
        vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
        vcf_readers.append(vcf_reader)
        vcf_reader.open()
        buffered_readers.append(BufferedReader(vcf_reader))

    return buffered_readers, vcf_readers

def _get_record_sample_data(vcf_record, format_tags):
    all_samples = {}
    for i in vcf_record.sample_dict:
        all_samples[i] = OrderedDict()

    for tag in format_tags:
        for i, sample_dict in vcf_record.sample_dict.items():
            if tag in sample_dict:
                all_samples[i][tag] = sample_dict[tag]
            else:
                all_samples[i][tag] = "."

    return all_samples

def _alter_record_fields(vcf_record, format_tags, coordinates):
    vcf_record.id = "."
    vcf_record.qual = "."
    vcf_record.filter = "."
    vcf_record.format_set = format_tags

    if len(coordinates) > 1:
        info = vcf_record.info.split(";") if vcf_record.info != "." else []
        info.append("JQ_MULT_ALT_LOCUS")
        vcf_record.info = ";".join(info)

    vcf_record.sample_dict = _get_record_sample_data(vcf_record, format_tags)

    return vcf_record

def add_subparser(subparser):
    #pylint: disable=C0301
    parser = subparser.add_parser("merge2", help="Accepts a directory of VCFs and returns a single merged VCF file.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to output variant-level VCF file")
    parser.add_argument("-a", "--allow_inconsistent_sample_sets", action="store_true", default=False, help="Allow inconsistent sample sets across callers. Not recommended.")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def execute(args, execution_context):
    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)

    all_metaheaders = []
    coordinate_set = set()
    input_files = sorted(glob.glob(os.path.join(input_path, "*.vcf")))

    file_writer = vcf.FileWriter(output_path)
    file_writer.open()

    count = 1
    samples = []
    for input_file in input_files:
        vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
        all_metaheaders = _produce_merged_metaheaders(vcf_reader,
                                                      all_metaheaders,
                                                      count)

        split_column_header = vcf_reader.column_header.split("\t")
        sample_cols = split_column_header[9:]
        samples.extend([vcf_reader.file_name + "|" +  i for i in sample_cols])
        column_header = split_column_header[0:9]

        coordinate_set = _add_to_coordinate_set(vcf_reader, coordinate_set)
        count += 1

    column_header.extend(samples)

    _write_metaheaders(file_writer,
                      all_metaheaders,
                      column_header,
                      execution_context)

    format_tags = _extract_format_ids(all_metaheaders)
    sorted_coordinates = _sort_coordinate_set(coordinate_set)

    buffered_readers, vcf_readers = _create_reader_lists(input_files)

    for coordinate in sorted_coordinates:
        for reader in buffered_readers:
            if reader.current_record:
                current_record = reader.current_record
                reader.check_current_value(coordinate)
                current_record = _alter_record_fields(current_record,
                                                      format_tags,
                                                      coordinate)
                line = current_record.asText()

        file_writer.write(line)

    for vcf_reader in vcf_readers:
        vcf_reader.close()

    file_writer.close()

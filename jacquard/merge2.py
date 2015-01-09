# pylint: disable=missing-docstring
from collections import defaultdict, OrderedDict
import glob
import os
import utils as utils
import vcf as vcf

# pylint: disable=too-few-public-methods
# This class must capture the state of the incoming iterator and provide
# modified behavior based on data in that iterator. A small class works ok, but
# suspect there may be a more pythonic way to curry iterator in a partial
# function. (cgates)
class GenericBufferedReader(object):
    '''Accepts an iterator and returns element (advancing current element) if
    requested element equals next element in collection and None otherwise.
    Never returns StopIteration so cannot be used as iteration control-flow.'''
    def __init__(self, iterator):
        self._iterator = iterator
        self._current_element = self._iterator.next()

    def get_if_equals(self, requested_element):
        result = None
        if requested_element == self._current_element:
            result = self._current_element
            self._current_element = self._get_next()
        return result

    def _get_next(self):
        try:
            return self._iterator.next()
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
    mult_alt_header = '##INFO=<ID=JQ_MULT_ALT_LOCUS,Number=0,Type=Flag,'\
                      'Description="dbSNP Membership",Source="Jacquard",'\
                      'Version="{}">'.format(utils.__version__)
    if mult_alt_header not in all_meta_headers:
        all_meta_headers.append(mult_alt_header)

    vcf_reader.close()

    return all_meta_headers

def _write_metaheaders(file_writer, all_headers, execution_context):
    column_header = all_headers.pop()
    all_headers.extend(execution_context)

    file_writer.write("\n".join(all_headers) + "\n")
    file_writer.write(column_header + "\n")

def _create_reader_lists(input_files):
    buffered_readers = []
    vcf_readers = []

    for input_file in input_files:
        vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
        vcf_readers.append(vcf_reader)
        vcf_reader.open()

        records = vcf_reader.vcf_records(qualified=True)
        buffered_readers.append(GenericBufferedReader(records))

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

def _build_coordinates(vcf_readers):
#     sample_list = []
    coordinate_set = OrderedDict()
    mult_alts = defaultdict(set)

    for vcf_reader in vcf_readers:
        try:
            vcf_reader.open()
#             sample_list.extend(vcf_reader.sample_names)
            for vcf_record in vcf_reader.vcf_records():
                coordinate_set[(vcf_record.get_empty_record())] = 0
                mult_alts[(vcf_record.chrom, vcf_record.pos, vcf_record.ref)]\
                    .add(vcf_record.alt)
        finally:
            vcf_reader.close()

    for vcf_record in coordinate_set:
        alts_for_this_locus = mult_alts[vcf_record.chrom,
                                        vcf_record.pos,
                                        vcf_record.ref]
        if len(alts_for_this_locus) > 1:
            #TODO: (cgates) move this logic inside VcfRecord
            info = vcf_record.info.split(";") if vcf_record.info != "." else []
            info.append("JQ_MULT_ALT_LOCUS")
            vcf_record.info = ";".join(info)

    coordinate_list = coordinate_set.keys()
    coordinate_list.sort()
#     sample_list.sort()

    return coordinate_list

def _build_merged_record(coordinate, vcf_records, all_sample_names):
    all_tags = set()
    sparse_matrix = {}

    for record in vcf_records:
        for sample, tags in record.sample_tag_values.items():
            sparse_matrix[sample] = {}
            for tag, value in tags.items():
                all_tags.add(tag)
                sparse_matrix[sample][tag] = value

    full_matrix = OrderedDict()
    for sample in all_sample_names:
        full_matrix[sample] = OrderedDict()
        for tag in sorted(all_tags):
            try:
                full_matrix[sample][tag] = sparse_matrix[sample][tag]
            except KeyError:
                full_matrix[sample][tag] = "."

    merged_record = vcf.VcfRecord(coordinate.chrom,
                                  coordinate.pos,
                                  coordinate.ref,
                                  coordinate.alt,
                                  coordinate.id,
                                  coordinate.qual,
                                  coordinate.filter,
                                  coordinate.info,
                                  sample_tag_values=full_matrix)

    return merged_record

def _pull_matching_records(coordinate, buffered_readers):
    vcf_records = []
    for reader in buffered_readers:
        record = reader.get_if_equals(coordinate)
        if record:
            vcf_records.append(record)

    return vcf_records

def _merge_records(coordinates, buffered_readers, writer, all_sample_names):
    for coordinate in coordinates:
        vcf_records = _pull_matching_records(coordinate, buffered_readers)
        merged_record = _build_merged_record(coordinate,
                                             vcf_records,
                                             all_sample_names)
        writer.write(merged_record.asText())

#TODO (cgates): Adjust to make better use of VcfReader;
#never parse outside of VcfReader/VcfRecord!
def _process_inputs(input_files):
    all_headers = []
    count = 1
    all_sample_names = []

    for input_file in input_files:
        vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
        all_headers = _produce_merged_metaheaders(vcf_reader,
                                                  all_headers,
                                                  count)
        column_header = vcf_reader.column_header.split("\t")[0:9]
        all_sample_names.extend([vcf_reader.file_name + "|" +  i for i in vcf_reader.sample_names])
        count += 1

    all_sample_names.sort()
    column_header.extend(all_sample_names)
    all_headers.append("\t".join(column_header))

    return all_headers, all_sample_names

def add_subparser(subparser):
    #pylint: disable=line-too-long
    parser = subparser.add_parser("merge2", help="Accepts a directory of VCFs and returns a single merged VCF file.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to output variant-level VCF file")
    parser.add_argument("-a", "--allow_inconsistent_sample_sets", action="store_true", default=False, help="Allow inconsistent sample sets across callers. Not recommended.")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def execute(args, execution_context):
    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)
    input_files = sorted(glob.glob(os.path.join(input_path, "*.vcf")))
    file_writer = vcf.FileWriter(output_path)
    file_writer.open()

    headers, all_sample_names = _process_inputs(input_files)

    _write_metaheaders(file_writer,
                       headers,
                       execution_context)
    buffered_readers, vcf_readers = _create_reader_lists(input_files)
    coordinates = _build_coordinates(vcf_readers)
    _merge_records(coordinates, buffered_readers, file_writer, all_sample_names)

    for vcf_reader in vcf_readers:
        vcf_reader.close()
    file_writer.close()

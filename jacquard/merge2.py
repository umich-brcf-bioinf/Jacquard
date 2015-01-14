# pylint: disable=missing-docstring, too-many-locals, too-few-public-methods
from __future__ import print_function, absolute_import
from collections import defaultdict, OrderedDict
import glob
import jacquard.utils as utils
import jacquard.vcf as vcf
import os
import re
import jacquard.logger as logger

MULT_ALT_HEADER = ('##INFO=<ID=JQ_MULT_ALT_LOCUS,Number=0,Type=Flag,'
                   'Description="dbSNP Membership",Source="Jacquard",'
                   'Version="{}">').format(utils.__version__)


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

def _merge_existing_metaheaders(vcf_readers, tags_to_keep):
    all_meta_headers = utils.OrderedSet()
    all_tags_to_keep = []
    visited_metaheaders = set()

    all_meta_headers.add(MULT_ALT_HEADER)

    for vcf_reader in vcf_readers:
        for tag_regex in tags_to_keep:
            for tag, format_metaheader in vcf_reader.format_metaheaders.items():
                visited_metaheaders.add(format_metaheader)
                if re.match(tag_regex+"$", tag):
                    all_tags_to_keep.append(tag)
                    all_meta_headers.add(format_metaheader)

#TODO: (jebene) - Only include JQ_MULT_ALT_LOCUS metaheaders if there is at least one mult-alt in the file
        for tag, info_metaheader in vcf_reader.info_metaheaders.items():
            visited_metaheaders.add(info_metaheader)
            if tag == "JQ_MULT_ALT_LOCUS":
                all_tags_to_keep.append(tag)
                all_meta_headers.add(info_metaheader)

#TODO: (jebene) - Since the FILTER fields will always be null, we could re-think the way we're doing this.
        for tag, filter_metaheader in vcf_reader.filter_metaheaders.items():
            visited_metaheaders.add(filter_metaheader)

        orig_metaheaders = set(vcf_reader.metaheaders)
        remaining_metaheaders = orig_metaheaders.difference(visited_metaheaders)
#TODO: (jebene) - add an update() method to utils.OrderedSet() to avoid this iteration
        for meta_header in remaining_metaheaders:
            all_meta_headers.add(meta_header)

    return all_meta_headers, all_tags_to_keep

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
    coordinate_set = OrderedDict()
    mult_alts = defaultdict(set)
    error = 0

    for vcf_reader in vcf_readers:
        previous_record = None
        try:
            vcf_reader.open()

            for vcf_record in vcf_reader.vcf_records():
                if previous_record and vcf_record < previous_record:
                    logger.error("VCF File [{}] is not sorted."
                                 .format(vcf_reader.file_name))
                    error = 1
                previous_record = vcf_record
                coordinate_set[(vcf_record.get_empty_record())] = 0
                mult_alts[(vcf_record.chrom, vcf_record.pos, vcf_record.ref)]\
                    .add(vcf_record.alt)
        finally:
            vcf_reader.close()

    if error:
        raise utils.JQException("One or more VCF files were not sorted. "
                                "Review inputs and try again.")

#     mult_alt_in_file = 0
    for vcf_record in coordinate_set:
        alts_for_this_locus = mult_alts[vcf_record.chrom,
                                        vcf_record.pos,
                                        vcf_record.ref]
        if len(alts_for_this_locus) > 1:
            vcf_record.add_info_field("JQ_MULT_ALT_LOCUS")
#             mult_alt_in_file = 1

    coordinate_list = coordinate_set.keys()
    coordinate_list.sort()

#     return coordinate_list, mult_alt_in_file
    return coordinate_list


def _build_merged_record(coordinate,
                         vcf_records,
                         all_sample_names,
                         tags_to_keep):

    all_tags = set()
    sparse_matrix = {}

    for record in vcf_records:
        for sample, tags in record.sample_tag_values.items():
            if sample not in sparse_matrix:
                sparse_matrix[sample] = {}
            for tag, value in tags.items():
                if tag in tags_to_keep:
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

def _merge_records(coordinates,
                   buffered_readers,
                   writer,
                   all_sample_names,
                   tags_to_keep):

    for coordinate in coordinates:
        vcf_records = _pull_matching_records(coordinate, buffered_readers)
        merged_record = _build_merged_record(coordinate,
                                             vcf_records,
                                             all_sample_names,
                                             tags_to_keep)
        writer.write(merged_record.asText())


def _process_headers(vcf_readers, tags_to_keep):
    all_sample_names = set()
    patient_to_file = defaultdict(list)

    (all_headers,
     all_tags_to_keep) = _merge_existing_metaheaders(vcf_readers, tags_to_keep)

    for vcf_reader in vcf_readers:
        patient = vcf_reader.file_name.split(".")[0]
        patient_samples = [patient + "|" +  i for i in vcf_reader.sample_names]

        for patient_sample in patient_samples:
            all_sample_names.add(patient_sample)
            patient_to_file[patient_sample].append(vcf_reader.file_name)

        #TODO: (cgates) Could we use a constant instead of parsing inside a loop?
        column_header = vcf_reader.column_header.split("\t")[0:9]

    patient_to_file_ordered = OrderedDict(sorted(patient_to_file.items(),
                                                 key=lambda t: t[0]))

    merge_metaheaders = _build_merge_metaheaders(patient_to_file_ordered)
    for metaheader in merge_metaheaders:
        all_headers.add(metaheader)

    sorted_all_sample_names = utils.NaturalSort(list(all_sample_names)).sorted
    column_header.extend(sorted_all_sample_names)
    all_headers.add("\t".join(column_header))

    return list(all_headers), sorted_all_sample_names, all_tags_to_keep


def _build_merge_metaheaders(patient_to_file):
    metaheaders = []
    metaheader_fmt = "##jacquard.merge.sample=<Column={0},Name={1},Source={2}>"
    for i, patient in enumerate(patient_to_file):
        filenames_string = "|".join(patient_to_file[patient])
        metaheader = metaheader_fmt.format(i+1, patient, filenames_string)
        metaheaders.append(metaheader)
    return metaheaders


def add_subparser(subparser):
    #pylint: disable=line-too-long
    parser = subparser.add_parser("merge2", help="Accepts a directory of VCFs and returns a single merged VCF file.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to output variant-level VCF file")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")
    parser.add_argument("--include_format_tags", dest='tags', help="Comma-separated list of regexs for format tags to include in output. Defaults to all JQ tags.")

def execute(args, execution_context):
    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)
    tags_to_keep = args.tags.split(",") if args.tags else ["JQ_.*"]

    input_files = sorted(glob.glob(os.path.join(input_path, "*.vcf")))
    file_writer = vcf.FileWriter(output_path)
    file_writer.open()

    buffered_readers, vcf_readers = _create_reader_lists(input_files)
    headers, all_sample_names, all_tags_to_keep = _process_headers(vcf_readers,
                                                                   tags_to_keep)

    _write_metaheaders(file_writer,
                       headers,
                       execution_context)

    coordinates = _build_coordinates(vcf_readers)

    _merge_records(coordinates,
                   buffered_readers,
                   file_writer,
                   all_sample_names,
                   all_tags_to_keep)

    for vcf_reader in vcf_readers:
        vcf_reader.close()
    file_writer.close()

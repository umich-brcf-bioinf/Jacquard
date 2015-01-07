# pylint: disable=missing-docstring
from collections import defaultdict, OrderedDict
import glob
import os
import re
from sets import Set

import utils as utils
import vcf as vcf

class BufferedReader(object):
    '''Accepts an iterator of ordered values, returns element if requested
    element matches and None otherwise. Never returns StopIteration so cannot
    be used as iteration control-flow.'''
    def __init__(self, reader):
        self.reader = reader
        self.base_iterator = reader.vcf_records()
        self.current_record = self.base_iterator.next()

    def extract(self, element):
        sample_dict = {}
        if element:
            for i, sample_name in enumerate(self.reader.samples):
                key = self.reader.file_name + "|" + sample_name
                sample_dict[key] = element.sample_dict[i]
        return sample_dict

    #TODO (cgates): Move this inside VcfRecord
    def get_sample_info(self, requested_element):
        result = self.extract(None)
        if requested_element == self.current_record:
            result = self.extract(self.current_record)
            self.current_record = self._get_next()
        return result

    def get_if_equals(self, requested_element):
        if requested_element == self.current_record:
            result = self.extract(self.current_record)
            self.current_record = self._get_next()
        return result

    def _get_next(self):
        try:
            return self.base_iterator.next()
        except StopIteration:
            return None

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

    def check_current_element(self):
        return self._current_element
    
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

#TODO (cgates): Remove this and associated tests
def _add_to_coordinate_set(vcf_reader, coordinate_set):
    vcf_reader.open()
    for vcf_record in vcf_reader.vcf_records():
        coordinate_set.add(vcf_record.get_empty_record())

    vcf_reader.close()

    return coordinate_set

#TODO (cgates): Remove this and associated tests
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

#pylint: disable=line-too-long
def add_subparser(subparser):
    parser = subparser.add_parser("merge2", help="Accepts a directory of VCFs and returns a single merged VCF file.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to output variant-level VCF file")
    parser.add_argument("-a", "--allow_inconsistent_sample_sets", action="store_true", default=False, help="Allow inconsistent sample sets across callers. Not recommended.")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def _build_coordinates(vcf_readers):
    coordinate_set = OrderedDict()
    mult_alts = defaultdict(set)

    for vcf_reader in vcf_readers:
        try:
            vcf_reader.open()
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
            #TODO (cgates): move this logic inside VcfRecord
            info = vcf_record.info.split(";") if vcf_record.info != "." else []
            info.append("JQ_MULT_ALT_LOCUS")
            vcf_record.info = ";".join(info)

    coordinate_list = coordinate_set.keys()
    coordinate_list.sort()
    return coordinate_list

def _build_tag_set(vcf_records):
    sample_tag_unordered = {}
    
    tag_set = set()
    for record in vcf_records:
        sample_tag_unordered.update(record.sample_tag_values)
        first_sample = record.sample_tag_values.keys()[0]
        for tag in record.sample_tag_values[first_sample]: 
            tag_set.add(tag)
    
    return sample_tag_unordered, sorted(list(tag_set))
    
def _sort_sample_tag_values(sample_tag_unordered, tag_tracker_list, all_sample_names):
    sample_tag_ordered = OrderedDict()
    
    for sample in sample_tag_unordered.keys():
        all_sample_names.add(sample)
        
    for sample in sorted(list(all_sample_names)):
        if sample in sample_tag_unordered:
            sorted_tags = OrderedDict()
            for tag in tag_tracker_list:
                sorted_tags[tag] = "."
                if tag in sample_tag_unordered[sample]:
                    tag_value = sample_tag_unordered[sample][tag]
                    sorted_tags[tag] = tag_value
                    
            sample_tag_ordered[sample] = sorted_tags
            
        else:
            sample_tag_ordered[sample] = OrderedDict()
            
    return sample_tag_ordered
        
def _build_merged_record(coordinate, vcf_records, all_sample_names):
    sample_tag_unordered, tag_tracker_list = _build_tag_set(vcf_records)
    sample_tag_ordered = _sort_sample_tag_values(sample_tag_unordered, tag_tracker_list, all_sample_names)
        
    merged_record = vcf.VcfRecord(coordinate.chrom, coordinate.pos, 
                                  coordinate.ref, coordinate.alt,
                                  coordinate.id, coordinate.qual, 
                                  coordinate.filter, coordinate.info,
                                  sample_tag_values = sample_tag_ordered)
    
    return merged_record 

# def _get_tag_sample_values(buffered_readers, merged_record):
#     tag_dict = defaultdict(dict)
#     for reader in buffered_readers:
#         record = reader.get_if_equals(merged_record)
#         if record:
#             sample_tag = record.sample_tag_values
#             for sample in sample_tag:
#                 for tag in sample_tag[sample]:
#                     value = sample_tag[sample][tag]
#                     tag_dict[tag].update({sample:value})
#             
#     return tag_dict

def _pull_matching_records(coordinate, buffered_readers):
    vcf_records = []
    for reader in buffered_readers:
        record = reader.get_if_equals(coordinate)
        if record:
            vcf_records.append(record)
            
    return vcf_records

def _merge_records(coordinates, buffered_readers, writer):
    all_sample_names = set()
    for reader in buffered_readers:
        record = reader.check_current_element()
        for sample in record.sample_tag_values:
            all_sample_names.add(sample)
        
    for coordinate in coordinates:
        vcf_records = _pull_matching_records(coordinate, buffered_readers)
        merged_record = _build_merged_record(coordinate, vcf_records, all_sample_names)
        writer.write(merged_record.asText())
        ##delete from coordinates

#         for reader in buffered_readers:
#             total_sample_dict.update(reader.get_sample_info(dest_record))
#         dest_record.set_sample_dict(total_sample_dict)
#         writer.write(dest_record.asText())

#TODO (cgates): Rewrite this to use build_coordinates
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
                reader.get_sample_info(coordinate)
                line = current_record.asText()

        file_writer.write(line)

    for vcf_reader in vcf_readers:
        vcf_reader.close()

    file_writer.close()

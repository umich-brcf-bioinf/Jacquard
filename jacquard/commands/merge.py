"""Merges a set of VCF files into a single VCF file.

* Merge assumes the incoming VCFs are aligned to the same set of contigs.
* Merge assumes incoming VCFs names follow this pattern:
    patientIdentifier.*.vcf
    Specifically, the first element of the VCF file name should be the patient
    name; for example you mihght have this:
        patientA.mutect.vcf, patientA.strelka.vcf, patientA.varscan.vcf,
        patientB.mutect.vcf, patientB.strelka.vcf, patientB.varscan.vcf,
        etc.
* Merge assumes the sample names (i.e. the VCF sample column headers, typically
    TUMOR and NORMAL) are consistent across the input VCFs. (The preceding
    Jacquard command "translate" ensures this is true.)
* Each incoming VCF record is joined with other VCF records that share the same
    "coordinate", where coordinate is the (chrom, pos, ref, and alt).
* Calls for the same patient/sample are joined into a single column. Each output
    column is the patient name (prefix of the file name) plus the sample name
    from the column header.
* The merged file will have as many records as the distinct set of
    (chrom, pos, ref, alt) across all input files.
* Each variant record will have the minimal set of incoming format tags for
    that variant (i.e. the list of format tags is specific to each record).
* Incoming QUAL and INFO fields are ignored.
* By default, merge only includes "Jacquard" FORMAT tags, but other tags can be
    included through and optional a command line arg.
"""
from __future__ import print_function, absolute_import, division

from collections import defaultdict, OrderedDict
import glob
import os
import re

import natsort

import jacquard.logger as logger
import jacquard.utils as utils
from jacquard.vcf import FileWriter
import jacquard.vcf as vcf


_DEFAULT_INCLUDED_FORMAT_TAGS = ["JQ_.*"]
_MULT_ALT_TAG = "JQ_MULT_ALT_LOCUS"
_JQ_SOMATIC_TAG = "HC_SOM"
#TODO: (jebene) is the description of this really accurate?
_MULT_ALT_HEADER = ('##INFO=<ID={},Number=0,Type=Flag,'
                    'Description="dbSNP Membership">').format(_MULT_ALT_TAG)
_FILE_FORMAT = ["##fileformat=VCFv4.1"]
_FILE_OUTPUT_SUFFIX = "merged"

class _BufferedReader(object):
    """A look-ahead reader that expedites merging.

    Accepts an iterator and returns element (advancing current element) if
    requested element equals next element in collection and None otherwise.
    Never returns StopIteration so cannot be used as iteration control-flow.

    This behavior is useful only because VcfRecord equality is based on their
    coordinate (chrom, pos, ref, alt), so by iterating over a list of
    coordinates, you can either park on your next record or return the
    current record and advance the reader to the next.

    Each incoming VcfReader is wrapped in a a BufferedReader and readers are
    advanced to the next reader when their current coordinate is requested.
    This approach avoids reading all VCFs into memory, but does require a file
    handle for each VCF you are merging and also requires that VcfReaders are
    sorted identically.

    Stylistic note:
    This class must capture the state of the incoming iterator and provide
    modified behavior based on data in that iterator. A small class works ok,
    but suspect there may be a more pythonic way to curry iterator in a
    partial function. Uncertain if that would be clearer/simpler. [cgates]
    """
    #pylint: disable=too-few-public-methods
    def __init__(self, iterator):
        self._iterator = iterator
        try:
            self._current_element = next(self._iterator)
        except StopIteration:
            self._current_element = self._iterator

    def next_if_equals(self, requested_element):
        result = None
        if requested_element == self._current_element:
            result = self._current_element
            self._current_element = self._get_next()
        return result

    def _get_next(self):
        try:
            return next(self._iterator)
        except StopIteration:
            return None

def _build_format_tags(format_tag_regex, vcf_readers):
    retained_tags = set()
    regexes_used = set()
    for vcf_reader in vcf_readers:
        for tag_regex in format_tag_regex:
            for tag in vcf_reader.format_metaheaders:
                if re.match(tag_regex + "$", tag):
                    retained_tags.add(tag)
                    regexes_used.add(tag_regex)

    if len(retained_tags) == 0:
        msg = ("The specified format tag regex [{}] would exclude all format "
               "tags. Review inputs/usage and try again")
        raise utils.JQException(msg, format_tag_regex)

    unused_regexes = set(format_tag_regex).difference(regexes_used)
    if unused_regexes:
        for unused_regex in unused_regexes:
            msg = ("In the specified list of regexes {}, the regex [{}] does "
                   "not match any format tags; this expression may be "
                   "irrelevant.")
            logger.warning(msg, format_tag_regex, unused_regex)

    return sorted(list(retained_tags))

def _compile_metaheaders(incoming_headers,
                         vcf_readers,
                         all_sample_names,
                         contigs_to_keep,
                         format_tags_to_keep,
                         info_tags_to_keep):
    #pylint: disable=too-many-arguments
    ordered_metaheaders = list(incoming_headers)
    all_info_metaheaders = {}
    all_format_metaheaders = {}
    all_contig_metaheaders = {}

    for vcf_reader in vcf_readers:
        all_info_metaheaders.update(vcf_reader.info_metaheaders)
        all_format_metaheaders.update(vcf_reader.format_metaheaders)
        all_contig_metaheaders.update(vcf_reader.contig_metaheaders)

    all_info_metaheaders[_MULT_ALT_TAG] = _MULT_ALT_HEADER

    if all_contig_metaheaders:
        for contig_id in contigs_to_keep:
            ordered_metaheaders.append(all_contig_metaheaders[contig_id])
    for tag in info_tags_to_keep:
        ordered_metaheaders.append(all_info_metaheaders[tag])
    for tag in format_tags_to_keep:
        ordered_metaheaders.append(all_format_metaheaders[tag])

    sorted_metaheaders = utils.sort_metaheaders(ordered_metaheaders)

    column_header = vcf_readers[0].column_header.split("\t")[0:9]
    column_header.extend(all_sample_names)
    sorted_metaheaders.append("\t".join(column_header))

    return sorted_metaheaders

def _write_metaheaders(file_writer, all_headers):
    file_writer.write("\n".join(all_headers) + "\n")

def _create_reader_lists(file_readers):
    buffered_readers = []
    vcf_readers = []

    for file_reader in file_readers:
        vcf_reader = vcf.VcfReader(file_reader)
        vcf_readers.append(vcf_reader)
        vcf_reader.open()

        records = vcf_reader.vcf_records(qualified=True)
        buffered_readers.append(_BufferedReader(records))

    return buffered_readers, vcf_readers

def _is_somatic_variant(record):
    for sample in record.sample_tag_values:
        for tag in record.sample_tag_values[sample]:
            is_somatic = record.sample_tag_values[sample][tag] == "1"
            if re.search(_JQ_SOMATIC_TAG, tag) and is_somatic:
                return 1

def _build_coordinates(vcf_readers, args):
    passed_variants = args.include_variants == "passed"
    any_passed_loci = args.include_loci == "any_passed"
    somatic_variants = args.include_variants == "somatic"
    any_somatic_loci = args.include_loci == "any_somatic"

    coordinate_set = set()
    mult_alts = defaultdict(set)

    for i, vcf_reader in enumerate(vcf_readers):
        logger.info("Reading [{}] ({}/{})",
                    os.path.basename(vcf_reader.file_name),
                    i + 1,
                    len(vcf_readers))
        try:
            vcf_reader.open()

            for vcf_record in vcf_reader.vcf_records():
                if passed_variants or any_passed_loci:
                    if vcf_record.filter == "PASS":
                        coordinate_set.add(vcf_record.get_empty_record())
                elif somatic_variants or any_somatic_loci:
                    if _is_somatic_variant(vcf_record):
                        coordinate_set.add(vcf_record.get_empty_record())
                else:
                    if "JQ_EXCLUDE" not in vcf_record.filter:
                        coordinate_set.add(vcf_record.get_empty_record())
                ref_alt = vcf_record.ref, vcf_record.alt
                locus = vcf_record.chrom, vcf_record.pos
                mult_alts[locus].add(ref_alt)
        finally:
            vcf_reader.close()

    for vcf_record in coordinate_set:
        ref_alts_for_this_locus = mult_alts[vcf_record.chrom, vcf_record.pos]
        inferred_mult_alt = len(ref_alts_for_this_locus) > 1
        explicit_mult_alt = "," in next(iter(ref_alts_for_this_locus))[1]
        if inferred_mult_alt or explicit_mult_alt:
            vcf_record.add_info_field(_MULT_ALT_TAG)

    return sorted(list(coordinate_set))

def _write_headers(reader, file_writer):
    headers = reader.metaheaders
    headers.append(reader.column_header)

    file_writer.write("\n".join(headers) + "\n")

def _sort_vcf(reader, temp_dir):
    logger.info("Sorting vcf [{}]", reader.file_name)
    vcf_records = []
    sorted_dir = os.path.join(temp_dir, "tmp")
    os.makedirs(sorted_dir)
    reader.open()
    for vcf_record in reader.vcf_records():
        vcf_records.append(vcf_record)

    reader.close()
    vcf_records.sort()
    writer = FileWriter(os.path.join(sorted_dir,
                                     reader.file_name))
    writer.open()
    writer.write("\n".join(reader.metaheaders) + "\n")
    writer.write(reader.column_header + "\n")
    for vcf_record in vcf_records:
        writer.write(vcf_record.text())

    writer.close()
    reader = vcf.VcfReader(vcf.FileReader(writer.output_filepath))
    return reader

def _get_unsorted_readers(vcf_readers):
    unsorted_readers = []
    for reader in vcf_readers:
        previous_record = None
        reader.open()
        for vcf_record in reader.vcf_records():
            if previous_record and vcf_record < previous_record:
                logger.debug("VCF file:chrom:pos [{}:{}:{}] is out of order"
                             .format(reader.file_name,
                                     vcf_record.chrom,
                                     vcf_record.pos))
                unsorted_readers.append(reader)
                break
            else:
                previous_record = vcf_record
        reader.close()
    return unsorted_readers

def _sort_readers(vcf_readers, temp_dir):
    unsorted_readers = _get_unsorted_readers(vcf_readers)
    sorted_readers = []
    for reader in vcf_readers:
        if reader in unsorted_readers:
            reader = _sort_vcf(reader, temp_dir)
        sorted_readers.append(reader)
    return sorted_readers


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
            for tag, value in list(tags.items()):
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
                                  coordinate.vcf_id,
                                  coordinate.qual,
                                  coordinate.filter,
                                  coordinate.info,
                                  sample_tag_values=full_matrix)

    return merged_record

#TODO: (jebene) alter logic to reduce number of branches
def _pull_matching_records(args, coordinate, buffered_readers):
    passed_variants = args.include_variants == "passed"
    all_passed_loci = args.include_loci == "all_passed"
    somatic_variants = args.include_variants == "somatic"
    all_somatic_loci = args.include_loci == "all_somatic"

    vcf_records = []
    for buffered_reader in buffered_readers:
        record = buffered_reader.next_if_equals(coordinate)
        if record:
            if passed_variants:
                if record.filter == "PASS":
                    vcf_records.append(record)
            elif all_passed_loci:
                if record.filter == "PASS":
                    vcf_records.append(record)
                else:
                    return []
            elif somatic_variants:
                if _is_somatic_variant(record):
                    vcf_records.append(record)
            elif all_somatic_loci:
                if _is_somatic_variant(record):
                    vcf_records.append(record)
                else:
                    return []
            else:
##accounts for any_passed/any_somatic because build_coordinates handled it
                vcf_records.append(record)

    return vcf_records

def _merge_records(args,
                   coordinate,
                   buffered_readers,
                   all_sample_names,
                   tags_to_keep):

    vcf_records = _pull_matching_records(args, coordinate, buffered_readers)
    if vcf_records:
        merged_record = _build_merged_record(coordinate,
                                             vcf_records,
                                             all_sample_names,
                                             tags_to_keep)
        return merged_record
    return 0

def _build_sample_list(vcf_readers):
    all_sample_names = set()
    patient_to_file = defaultdict(list)

    for vcf_reader in vcf_readers:
        for sample_name in vcf_reader.qualified_sample_names:
            all_sample_names.add(sample_name)
            patient_to_file[sample_name].append(vcf_reader.file_name)
    sorted_sample_names = natsort.natsorted(all_sample_names)

    sorted_patient_to_file = OrderedDict()
    for sample_name in sorted_sample_names:
        sorted_patient_to_file[sample_name] = patient_to_file[sample_name]

    merge_metaheaders = _build_merge_metaheaders(sorted_patient_to_file)

    return sorted_sample_names, merge_metaheaders

def _build_info_tags(coordinates):
    all_info_tags = set()
    for record in coordinates:
        all_info_tags.update(record.info_dict.keys())
    ordered_tags = sorted(all_info_tags)

    return ordered_tags

def _build_contigs(coordinates):
    all_contigs = set()
    for record in coordinates:
        all_contigs.add(record.chrom)
    ordered_contigs = natsort.natsorted(all_contigs)

    return ordered_contigs

def _build_merge_metaheaders(patient_to_file):
    metaheaders = []
    metaheader_fmt = "##jacquard.merge.sample=<Column={0},Name={1},Source={2}>"
    for i, patient in enumerate(patient_to_file):
        filenames_string = "|".join(patient_to_file[patient])
        metaheader = metaheader_fmt.format(i + 1, patient, filenames_string)
        metaheaders.append(metaheader)

    return metaheaders

def _build_writers_to_readers(vcf_readers, output_path):
    writers_to_readers = {}
    basename, extension = os.path.splitext(output_path)
    new_filename = basename + ".merged" + extension

    file_writer = vcf.FileWriter(new_filename)
    writers_to_readers[file_writer] = vcf_readers

    return writers_to_readers

def _get_readers_per_patient(file_readers):
    readers_per_patient = defaultdict(list)
    for file_reader in file_readers:
        file_reader.open()
        patient = file_reader.file_name.split(".")[0]
        for line in file_reader.read_lines():
            caller_meta_header = "##jacquard.translate.caller"
            if line.startswith(caller_meta_header):
                readers_per_patient[patient].append(line.split("=")[1])
        file_reader.close()

    return OrderedDict(sorted(readers_per_patient.items()))

def _validate_consistent_samples(file_readers):
    readers_per_patient = _get_readers_per_patient(file_readers)

    all_callers = set()
    for callers in readers_per_patient.values():
        all_callers.update(callers)
    warning = 0
    for patient, callers in readers_per_patient.items():
        missing_callers = set(all_callers).difference(set(callers))
        if missing_callers:
            warning = 1
            msg = "Sample [{}] is missing VCF(s): {}"
            logger.warning(msg,
                           patient,
                           sorted(list(missing_callers)))
    if warning:
        msg = "Some samples appear to be missing VCF(s)"
        logger.warning(msg)

def add_subparser(subparser):
    #pylint: disable=line-too-long
    parser = subparser.add_parser("merge", help="Accepts a directory of VCFs and returns a single merged VCF file.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to output variant-level VCF file")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")
    parser.add_argument("--include_format_tags", dest='tags', help="Comma-separated list of regexs for format tags to include in output. Defaults to all JQ tags.")
    parser.add_argument("--log_file", help="Log file destination")
    parser.add_argument("--include_variants", choices=["passed", "somatic"], help=("passed: Only include variants which passed their respective filter\n"
                                                                                   "somatic: Only include somatic variants"))

    parser.add_argument("--include_loci", choices=["any_passed", "all_passed", "any_soamtic", "all_somatic"], help=("any_passed: Include all variants at loci where at least one variant passed\n"
    # pylint: disable=line-too-long
                                                                                                                        "all_passed: Include all variants at loci where all variants passed\n"
                                                                                                                        "any_somatic: Include all variants at loci where at least one variant was somatic\n"
                                                                                                                        "all_somatic: Include all variants at loci where all variants were somatic"))

def _predict_output(args):
    desired_output_files = set([os.path.basename(args.output)])

    return desired_output_files

def report_prediction(args):
    return _predict_output(args)

def get_required_input_output_types():
    return ("directory", "file")

#TODO (cgates): Validate should actually validate
def validate_args(dummy):
    pass

def execute(args, execution_context):
    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)
    if args.tags:
        format_tag_regex = args.tags.split(",")
    else:
        format_tag_regex = _DEFAULT_INCLUDED_FORMAT_TAGS

    input_files = sorted(glob.glob(os.path.join(input_path, "*.vcf")))
    file_readers = [vcf.FileReader(i) for i in input_files]
    _validate_consistent_samples(file_readers)

    try:
        file_writer = vcf.FileWriter(output_path)
        file_writer.open()

        buffered_readers, vcf_readers = _create_reader_lists(file_readers)
        vcf_readers = _sort_readers(vcf_readers, output_path)
        all_sample_names, merge_metaheaders = _build_sample_list(vcf_readers)
        coordinates = _build_coordinates(vcf_readers, args)
        format_tags_to_keep = _build_format_tags(format_tag_regex, vcf_readers)
        info_tags_to_keep = _build_info_tags(coordinates)
        contigs_to_keep = _build_contigs(coordinates)
        incoming_headers = _FILE_FORMAT + execution_context + merge_metaheaders
        headers = _compile_metaheaders(incoming_headers,
                                       vcf_readers,
                                       all_sample_names,
                                       contigs_to_keep,
                                       format_tags_to_keep,
                                       info_tags_to_keep)

        _write_metaheaders(file_writer, headers)

        for coordinate in coordinates:
            merged_record = _merge_records(args,
                                           coordinate,
                                           buffered_readers,
                                           all_sample_names,
                                           format_tags_to_keep)
            if merged_record:
                file_writer.write(merged_record.text())

    finally:
        for vcf_reader in vcf_readers:
            vcf_reader.close()
        file_writer.close()

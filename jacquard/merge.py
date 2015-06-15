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
* By default, the merged file will have as many records as the distinct set of
    (chrom, pos, ref, alt) across all input files.
    Running Merge with "--include_variants" or "--include_loci" will reduce the
    number of records in the merged file to include only those which conform
    with the desired flag.
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
import textwrap

import jacquard.utils.logger as logger
import jacquard.utils.utils as utils
from jacquard.utils.vcf import FileWriter
import jacquard.utils.vcf as vcf


_DEFAULT_INCLUDED_FORMAT_TAGS = ["JQ_.*"]
_MULT_ALT_TAG = "JQ_MULT_ALT_LOCUS"
_JQ_SOMATIC_TAG = "HC_SOM"

_MULT_ALT_HEADER = ('##INFO=<ID={},Number=0,Type=Flag,'
                    'Description="More than one alt allele was seen at this '
                    'locus.">').format(_MULT_ALT_TAG)
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

class MergeVcfReader(vcf.VcfReader):
    def __init__(self, file_reader):
        super(self.__class__,self).__init__(file_reader)
        self.format_tags = {}

    def modify_metaheader(self, original_metaheader, transformed_tag):
        updated_metaheader = re.sub(r'(^##FORMAT=.*?[<,]ID=)([^,>]*)',
                                    r'\g<1>%s' % transformed_tag,
                                    original_metaheader)

        self.metaheaders.append(updated_metaheader)
        if original_metaheader in self.metaheaders:
            self.metaheaders.remove(original_metaheader)

    def store_format_tags(self, original_tag, new_tag):
        self.format_tags[original_tag] = new_tag

    @staticmethod
    def modify_format_tag(vcf_record, format_tags):
        for tags in list(vcf_record.sample_tag_values.values()):
            for original_tag, new_tag in list(format_tags.items()):
                if new_tag not in tags and original_tag in tags:
                    tags[new_tag] = tags[original_tag]
                    del tags[original_tag]

        return vcf_record

    def vcf_records(self, format_tags=None, qualified=False):
        """Generates parsed VcfRecord objects.

        Typically called in a for loop to process each vcf record in a
        VcfReader. VcfReader must be opened in advanced and closed when
        complete. Skips all headers.

        Args:
            qualified: When True, sample names are prefixed with file name

        Returns:
            Parsed VcfRecord

        Raises:
            StopIteration: when reader is exhausted.
            TypeError: if reader is closed.
        """
        if qualified:
            sample_names = self.qualified_sample_names
        else:
            sample_names = self.sample_names

        for line in self._file_reader.read_lines():
            if line.startswith("#"):
                continue
            vcf_record = vcf.VcfRecord.parse_record(line, sample_names)
            if format_tags:
                vcf_record = self.modify_format_tag(vcf_record, format_tags)
            yield vcf_record


class _Filter(object):
    """Interprets command-line args to initialize variant and locus filters
    """
    #pylint: disable=too-few-public-methods
    def __init__(self, args):
        self._cell_filters = {"all" : _Filter._include_cell_if_all,
                              "valid": _Filter._include_cell_if_valid,
                              "passed" : _Filter._include_cell_if_passed,
                              "somatic" : _Filter._include_cell_if_somatic}
        self._row_filters = {"all": _Filter._include_row_if_all,
                             "all_passed": _Filter._include_row_if_all_passed,
                             "at_least_one_passed": _Filter.\
                                                    _include_row_if_any_passed,
                             "all_somatic": _Filter._include_row_if_all_somatic,
                             "at_least_one_somatic": _Filter.\
                                                    _include_row_if_any_somatic}
        self._args = args
#         if args.include_all and (args.include_cells or args.include_rows):
#             msg = ('Unable to process command-line arguments. Neither '
#                    '--include_cells nor --include_rows can be specified if '
#                    '--include_all is specified.')
#             raise utils.UsageError(msg)

        if args.include_all:
            cell_filter_strategy = self._cell_filters["all"]
            row_filter_strategy = self._row_filters["all"]
        else:
            cell_filter_strategy = self._cell_filters[args.include_cells]
            row_filter_strategy = self._row_filters[args.include_rows]

        self._cell_filter_strategy = cell_filter_strategy
        self._row_filter_strategy = row_filter_strategy
        self.row_count = 0
        self.rows_excluded = 0
        self.cell_count = 0
        self.cells_excluded = 0

        self.excluded_breakdown = defaultdict(int)

    def include_cell(self, vcf_record):
        self.cell_count += 1
        include = self._cell_filter_strategy(self, vcf_record)
        if not include:
            self.cells_excluded += 1
        return include

    def include_row(self, vcf_records):
        self.row_count += 1
        include = self._row_filter_strategy(vcf_records)
        if not include:
            self.rows_excluded += 1
        return include

    def log_statistics(self):
        cell_msg = "{}% ({}) cells were excluded because (--include_cells={})"
        logger.info(cell_msg,
                    int(round(100 * self.cells_excluded / self.cell_count, 0)),
                    self.cells_excluded,
                    self._args.include_cells)

        row_msg = "{}% ({}) rows were excluded because (--include_rows={})"
        logger.info(row_msg,
                    int(round(100 * self.rows_excluded / self.row_count, 0)),
                    self.rows_excluded,
                    self._args.include_rows)

        for key, count in self.excluded_breakdown.items():
            msg = "{} cells were excluded with [{}]"
            logger.debug(msg, count, key)

    #TODO: (cgates/jebene): Should this be part of VcfRecord?
    @staticmethod
    def _is_somatic(record):
        for sample in record.sample_tag_values:
            for tag in record.sample_tag_values[sample]:
                is_somatic = record.sample_tag_values[sample][tag] == "1"
                if re.search(_JQ_SOMATIC_TAG, tag) and is_somatic:
                    return True
        return False

    def _include_cell_if_all(self, dummy):
        #pylint: disable=no-self-use
        return True

    def _include_cell_if_somatic(self, record):
        include = _Filter._is_somatic(record)
        if not include:
            self.excluded_breakdown["not somatic"] += 1
        return include

    def _include_cell_if_valid(self, record):
        include ="JQ_EXCLUDE" not in record.filter
        if not include:
            self.excluded_breakdown[record.filter] += 1
        return include

    def _include_cell_if_passed(self, record):
        include = record.filter == "PASS"
        if not include:
            self.excluded_breakdown[record.filter] += 1
        return include

    @staticmethod
    def _include_row_if_all(dummy):
        return True

    @staticmethod
    def _include_row_if_all_passed(records):
        for record in records:
            if record.filter != "PASS":
                return False
        return True

    @staticmethod
    def _include_row_if_any_passed(records):
        for record in records:
            if record.filter == "PASS":
                return True
        return False

    @staticmethod
    def _include_row_if_all_somatic(records):
        for record in records:
            if not _Filter._is_somatic(record):
                return False
        return True

    @staticmethod
    def _include_row_if_any_somatic(records):
        for record in records:
            if _Filter._is_somatic(record):
                return True
        return False

def _get_caller_name(merge_vcf_reader):
    for metaheader in merge_vcf_reader.metaheaders:
        tag = re.match("##jacquard.translate.caller=(.*)", metaheader)
        if tag:
            return tag.group(1)
    return False

def _alter_description(metaheader, caller_name):
    return re.sub(r'^##FORMAT=.*?[<,]Description="([,>]*)',
                  r'\g<0>[%s]: ' % caller_name,
                  metaheader)

def _get_format_tags(merge_vcf_readers):
    format_tags = defaultdict(list)
    caller_names = defaultdict(list)
    for merge_vcf_reader in merge_vcf_readers:
        caller_name = _get_caller_name(merge_vcf_reader)
        for tag, metahdr in list(merge_vcf_reader.format_metaheaders.items()):
            format_tags[tag].append(metahdr)
            caller_names[tag].append(caller_name)

    for tag, callers in caller_names.items():
        if len(set(callers)) != 1:
            for i, caller in enumerate(callers):
                format_tags[tag][i] = _alter_description(format_tags[tag][i],
                                                         caller)
        else:
            format_tags[tag] = [format_tags[tag][0]]

    return format_tags

def _disambiguate_format_tags(merge_vcf_readers, format_tags):
    for tag, metaheaders in list(format_tags.items()):
        for merge_vcf_reader in merge_vcf_readers:
            for i, metahdr in enumerate(list(metaheaders)):
                vcf_metaheaders = merge_vcf_reader.format_metaheaders.values()
                if metahdr in list(vcf_metaheaders):
                    if len(metaheaders) > 1:
                        new_tag = "JX{}_{}".format(i+1, tag)
                        format_tags[new_tag] = metaheaders

                        merge_vcf_reader.modify_metaheader(metahdr, new_tag)
                        merge_vcf_reader.store_format_tags(tag, new_tag)
                    else:
                        merge_vcf_reader.store_format_tags(tag, tag)

    return format_tags

def _build_format_tags(format_tag_regex, vcf_readers):
    retained_tags = set()
    regexes_used = set()

    for vcf_reader in vcf_readers:
        for tag_regex in format_tag_regex:
            for original_tag, new_tag in list(vcf_reader.format_tags.items()):
                if re.match(tag_regex + "$", original_tag):
                    retained_tags.add(new_tag)
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

def _create_vcf_readers(file_readers):
    merge_vcf_readers = []
    for file_reader in file_readers:
#         vcf_reader = vcf.VcfReader(file_reader)
        merge_vcf_reader = MergeVcfReader(file_reader)
        merge_vcf_readers.append(merge_vcf_reader)

    return merge_vcf_readers

def _create_buffered_readers(vcf_readers):
    buffered_readers = []
    for vcf_reader in vcf_readers:
        vcf_reader.open()
        records = vcf_reader.vcf_records(vcf_reader.format_tags, qualified=True)
        buffered_readers.append(_BufferedReader(records))

    return buffered_readers

def _build_coordinates(vcf_readers):
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

def _sort_vcf(reader, sorted_dir):
    vcf_records = []
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
    reader = MergeVcfReader(vcf.FileReader(writer.output_filepath))
    return reader

def _get_unsorted_readers(vcf_readers):
    unsorted_readers = []
    for i, reader in enumerate(vcf_readers):
        logger.info("Checking sort order of [{}] ({}/{})",
                    reader.file_name,
                    i+1,
                    len(vcf_readers)
                    )
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

def _sort_readers(vcf_readers, output_path):
    unsorted_readers = _get_unsorted_readers(vcf_readers)
    sorted_readers = []
    unsorted_count = 0
    if unsorted_readers:
        sorted_dir = os.path.join(os.path.dirname(output_path), "tmp")
        if not os.path.isdir(sorted_dir):
            os.makedirs(sorted_dir)

    for reader in vcf_readers:
        if reader in unsorted_readers:
            unsorted_count += 1
            logger.info("Sorting vcf [{}] ({}/{})",
                        reader.file_name,
                        unsorted_count,
                        len(unsorted_readers))
            reader = _sort_vcf(reader, sorted_dir)
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

def _pull_matching_records(filter_strategy, coordinate, buffered_readers):
    vcf_records = []
    for buffered_reader in buffered_readers:
        record = buffered_reader.next_if_equals(coordinate)
        if record and filter_strategy.include_cell(record):
            vcf_records.append(record)
    if filter_strategy.include_row(vcf_records):
        return vcf_records
    else:
        return []

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

def _validate_consistent_input(vcf_readers, include_all):
    if not include_all:
        untranslated_readers = []
        translated = False
        untranslated = False

        for vcf_reader in vcf_readers:
            if '##jacquard.translate.caller' in vcf_reader.metaheaders:
                translated = True
            else:
                untranslated = True
                untranslated_readers.append(vcf_reader)

        if translated and untranslated:
            msg = ("Some input VCFs [{}] were not translated by Jacquard. "
                    "Review input and/or use --included_all flag")\
                   .format(untranslated_readers)
            raise utils.UsageError(msg)


def _merge_records(vcf_readers,
                   coordinates,
                   filter_strategy,
                   all_sample_names,
                   format_tags_to_keep,
                   file_writer):
    #pylint: disable=too-many-arguments
    row_count = 0
    next_breakpoint = 0
    buffered_readers = _create_buffered_readers(vcf_readers)

    for coordinate in coordinates:
        vcf_records = _pull_matching_records(filter_strategy,
                                             coordinate,
                                             buffered_readers)

        if vcf_records:
            merged_record = _build_merged_record(coordinate,
                                                 vcf_records,
                                                 all_sample_names,
                                                 format_tags_to_keep)
            file_writer.write(merged_record.text())

        row_count += 1
        progress = 100 * row_count / len(coordinates)
        if progress > next_breakpoint:
            logger.info("Merging: {} rows processed ({}%)",
                        row_count,
                        next_breakpoint)
            next_breakpoint = 10 * int(progress/10) + 10

    logger.info("Merge complete: {} rows processed (100%)", row_count)
    filter_strategy.log_statistics()

def _get_format_tag_regex(args):
    if args.tags and args.include_all:
        msg =("Unable to process command-line arguments. --include_format_tags cannot "
              "be specified if --include_all is specified.")
        raise utils.UsageError(msg)

    if args.include_all:
        format_tag_regex = ['.*']
    elif args.tags:
        format_tag_regex = args.tags.split(",")
    else:
        format_tag_regex = _DEFAULT_INCLUDED_FORMAT_TAGS

    return format_tag_regex

def add_subparser(subparser):
    #pylint: disable=line-too-long

    parser = subparser.add_parser("merge",
                                  formatter_class=utils._JacquardHelpFormatter,
                                  usage=[" [--include_format_tags=JQ_.*]\n"
                                         "\t\t\t\t\t[--include_cells=valid]\n"
                                         "\t\t\t\t\t[--include_rows=at_least_one_somatic]"],
                                  help="Accepts a directory of VCFs and returns a single merged VCF file.")
    parser.add_argument("input", help="Directory containing VCF files. Other file types ignored")
    parser.add_argument("output", help="VCF file")
    parser.add_argument("--include_all", action="store_true", help=textwrap.fill("Useful when merging untranslated VCFs which have already been filtered to passing variants. Equivalent to --include_format_tags='.*' --include_cells=all --include_rows=all"))
    parser.add_argument("--include_cells",
                        choices=["all", "valid", "passed", "somatic"],
                        default="valid",
                        help=("all: Include all variants\n"
                              "valid: Only include valid variants (invalid variants become '.')\n"
                              "passed: Only include variants which passed their respective filter (failed variants become '.')\n"
                              "somatic: Only include somatic variants (non-somatic variants become '.')"),
                        metavar="")
    parser.add_argument("--include_rows",
                        choices=["all", "at_least_one_passed", "all_passed", "at_least_one_somatic", "all_somatic"],
                        default="at_least_one_somatic",
                        help=("all: Include all loci\n"
                              "at_least_one_passed: Include loci where at least one variant passed\n"
                              "at_least_one_somatic: Include loci where at least one variant was somatic\n"
                              "all_passed: Include loci where all variants passed\n"
                              "all_somatic: Include loci where all variants were somatic\n"),
                        metavar="")
    parser.add_argument("--include_format_tags", dest='tags', help=textwrap.fill("Comma-separated user-defined list of regular expressions for format tags to be included in output"), metavar="")
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")
    parser.add_argument("--log_file", help="Log file destination", metavar="")
    parser.add_argument("-v", "--verbose", action='store_true')

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
    filter_strategy = _Filter(args)
    format_tag_regex = _get_format_tag_regex(args)

    input_files = sorted(glob.glob(os.path.join(input_path, "*.vcf")))
    file_readers = [vcf.FileReader(i) for i in input_files]
    _validate_consistent_samples(file_readers)

    try:
        file_writer = vcf.FileWriter(output_path)
        file_writer.open()

#TODO: jebene: _build_format_tags, _built_sample_list, _build_info_tags,
#and _build_contigs behave differently. It seems like we could make the
#signatures of these methods more similar or even combine some methods to
#reduce excess iterations over the coordinates/vcf_readers
        merge_vcf_readers = _create_vcf_readers(file_readers)
        _validate_consistent_input(merge_vcf_readers, args.include_all)
        merge_vcf_readers = _sort_readers(merge_vcf_readers, output_path)
        format_tags = _get_format_tags(merge_vcf_readers)
        _disambiguate_format_tags(merge_vcf_readers, format_tags)
        format_tags_to_keep = _build_format_tags(format_tag_regex,
                                                 merge_vcf_readers)
        (all_sample_names,
         merge_metaheaders) = _build_sample_list(merge_vcf_readers)

        coordinates = _build_coordinates(merge_vcf_readers)
        info_tags_to_keep = _build_info_tags(coordinates)
        contigs_to_keep = _build_contigs(coordinates)
        incoming_headers = _FILE_FORMAT + execution_context + merge_metaheaders
        headers = _compile_metaheaders(incoming_headers,
                                       merge_vcf_readers,
                                       all_sample_names,
                                       contigs_to_keep,
                                       format_tags_to_keep,
                                       info_tags_to_keep)

        _write_metaheaders(file_writer, headers)

        _merge_records(merge_vcf_readers,
                       coordinates,
                       filter_strategy,
                       all_sample_names,
                       format_tags_to_keep,
                       file_writer)
    finally:
        for vcf_reader in merge_vcf_readers:
            vcf_reader.close()
        file_writer.close()

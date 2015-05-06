"""Adds summary tags/fields to a merged VCF file.

Collaborates with two summary "callers" to add INFO and FORMAT tags to each
variant record based on the presence of previously translated tags.
"""
from __future__ import print_function, absolute_import, division

import argparse
import os

import jacquard.utils.logger as logger
import jacquard.utils.utils as utils
import jacquard.utils.summarize_rollup_transform as summarize_caller
import jacquard.utils.summarize_zscore_transform as zscore_caller
import jacquard.utils.vcf as vcf


def _write_metaheaders(caller,
                       vcf_reader,
                       file_writer,
                       execution_context=None,
                       new_meta_headers=None):

    new_headers = list(vcf_reader.metaheaders)

    if execution_context:
        new_headers.extend(execution_context)
        new_headers.extend(caller.get_metaheaders())
    if new_meta_headers:
        new_headers.append(new_meta_headers)

    sorted_metaheaders = utils.sort_metaheaders(new_headers)
    sorted_metaheaders.append(vcf_reader.column_header)

    file_writer.write("\n".join(sorted_metaheaders) +"\n")

def _write_to_tmp_file(caller, vcf_reader, tmp_writer):
    vcf_reader.open()
    tmp_writer.open()

    try:
        _write_metaheaders(caller, vcf_reader, tmp_writer)
        logger.info("Adding summary tags for [{}]", vcf_reader.file_name)
        _add_tags(caller, vcf_reader, tmp_writer)

    finally:
        vcf_reader.close()
        tmp_writer.close()


def _write_zscores(caller,
                   metaheaders,
                   vcf_reader,
                   file_writer):

#TODO: (jebene) make zscores and tmp file follow the same pattern when writing
    try:
        file_writer.open()
        headers = list(metaheaders)
        headers.extend(vcf_reader.metaheaders)
        headers.extend(caller.metaheaders)
        sorted_metaheaders = utils.sort_metaheaders(headers)
        sorted_metaheaders.append(vcf_reader.column_header)
        file_writer.write("\n".join(sorted_metaheaders) +"\n")

        vcf_reader.open()
        for vcf_record in vcf_reader.vcf_records():
            line = caller.add_tags(vcf_record)
            file_writer.write(line)
    finally:
        vcf_reader.close()
        file_writer.close()

def _add_tags(caller, vcf_reader, file_writer):
    for vcf_record in vcf_reader.vcf_records():
        caller.add_tags(vcf_record)
        file_writer.write(vcf_record.text())

def add_subparser(subparser):
    # pylint: disable=line-too-long
    parser = subparser.add_parser("summarize", formatter_class=argparse.RawTextHelpFormatter, help="Accepts a Jacquard-merged VCF file and creates a new file, adding summary fields.")
    parser.add_argument("input", help="Path to Jacquard-merged VCF (or any VCF with Jacquard tags; e.g. JQ_SOM_MT)")
    parser.add_argument("output", help="Path to output VCf")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")
    parser.add_argument("--log_file", help="Log file destination")

def report_prediction(args):
    return set([os.path.basename(args.output)])

def get_required_input_output_types():
    return ("file", "file")

#TODO (cgates): Validate should actually validate
def validate_args(dummy):
    pass

def execute(args, execution_context):
    input_file = os.path.abspath(args.input)
    output = os.path.abspath(args.output)

    summary_caller = summarize_caller.SummarizeCaller()

    vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
    tmp_output_file = output + ".tmp"
    tmp_writer = vcf.FileWriter(tmp_output_file)

    _write_to_tmp_file(summary_caller, vcf_reader, tmp_writer)

    tmp_reader = vcf.VcfReader(vcf.FileReader(tmp_output_file))
    file_writer = vcf.FileWriter(output)

    logger.info("Calculating zscores")
    caller = zscore_caller.ZScoreCaller(tmp_reader)
    metaheaders = execution_context + summary_caller.get_metaheaders()
    _write_zscores(caller, metaheaders, tmp_reader, file_writer)

    os.remove(tmp_output_file)

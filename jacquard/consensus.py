from __future__ import print_function, absolute_import

import os

import jacquard.logger as logger
import jacquard.utils as utils
import jacquard.variant_callers.jacquard_consensus_caller as consensus_caller
import jacquard.variant_callers.jacquard_zscore_caller as zscore_caller
import jacquard.vcf as vcf


def _write_metaheaders(cons_helper,
                       vcf_reader,
                       file_writer,
                       execution_context=None,
                       new_meta_headers=None):

    new_headers = list(vcf_reader.metaheaders)

    if execution_context:
        new_headers.extend(execution_context)
        new_headers.extend(cons_helper.get_consensus_metaheaders())
    if new_meta_headers:
        new_headers.append(new_meta_headers)

    new_headers.append(vcf_reader.column_header)
    file_writer.write("\n".join(new_headers) +"\n")

def _write_to_tmp_file(cons_helper, vcf_reader, tmp_writer):
    vcf_reader.open()
    tmp_writer.open()

    try:
        _write_metaheaders(cons_helper, vcf_reader, tmp_writer)
        logger.info("Adding consensus tags for [{}]", vcf_reader.file_name)
        _add_consensus_tags(cons_helper, vcf_reader, tmp_writer)

    finally:
        vcf_reader.close()
        tmp_writer.close()


def _write_zscores(caller,
                   metaheaders,
                   vcf_reader,
                   file_writer):

    try:
        file_writer.open()
        headers = list(metaheaders)
        headers.extend(vcf_reader.metaheaders)
        headers.extend(caller.metaheaders)
        headers.append(vcf_reader.column_header)
        file_writer.write("\n".join(headers) +"\n")

        vcf_reader.open()
        for vcf_record in vcf_reader.vcf_records():
            line = caller.add_tags(vcf_record)
            file_writer.write(line)
    finally:
        vcf_reader.close()
        file_writer.close()

def _add_consensus_tags(cons_helper, vcf_reader, file_writer):
    for vcf_record in vcf_reader.vcf_records():
        file_writer.write(cons_helper.add_tags(vcf_record))

def add_subparser(subparser):
    # pylint: disable=line-too-long
    parser = subparser.add_parser("consensus", help="Accepts a Jacquard-merged VCf file and creates a new file, adding consensus fields.")
    parser.add_argument("input", help="Path to Jacquard-merged VCF (or any VCF with Jacquard tags (e.g. JQ_SOM_MT)")
    parser.add_argument("output", help="Path to output VCf")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def _predict_output(args):
    return set([os.path.basename(args.output)])

def report_prediction(args):
    return _predict_output(args)

def get_required_input_output_types():
    return ("file", "file")

def execute(args, execution_context):
    input_file = os.path.abspath(args.input)
    output = os.path.abspath(args.output)

    extension = os.path.splitext(os.path.basename(input_file))[1]
    if not os.path.isfile(input_file) or extension != ".vcf":
        raise utils.JQException("Input file [{}] must be a VCF file.",
                                input_file)

    extension = os.path.splitext(os.path.basename(output))[1]
    if extension:
        if extension != ".vcf":
            raise utils.JQException("Output file [{}] must be a VCF file.",
                                    output)
        else:
            utils.validate_directories(os.path.dirname(input_file),
                                       os.path.dirname(output))
            output_file = output
    else:
        utils.validate_directories(os.path.dirname(input_file), output)
        output_file = os.path.join(output, "consensus.vcf")

    cons_helper = consensus_caller.ConsensusCaller()

    vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
    tmp_output_file = output_file + ".tmp"
    tmp_writer = vcf.FileWriter(tmp_output_file)

    _write_to_tmp_file(cons_helper, vcf_reader, tmp_writer)

    tmp_reader = vcf.VcfReader(vcf.FileReader(tmp_output_file))
    file_writer = vcf.FileWriter(output_file)

    logger.info("Calculating zscores")
    caller = zscore_caller.ZScoreCaller(tmp_reader)
    metaheaders = execution_context + cons_helper.get_consensus_metaheaders()
    _write_zscores(caller, metaheaders, tmp_reader, file_writer)

    os.remove(tmp_output_file)

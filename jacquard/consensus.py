from __future__ import print_function, absolute_import
import os

from jacquard.variant_callers import consensus_helper
import jacquard.vcf as vcf
import jacquard.utils as utils
import jacquard.logger as logger

def _write_metaheaders(cons_helper,
                       vcf_reader,
                       file_writer,
                       execution_context=0,
                       new_meta_headers=0):

    new_headers = vcf_reader.metaheaders

    if execution_context:
        new_headers.extend(execution_context)
        new_headers.extend(cons_helper.get_consensus_metaheaders())
    if new_meta_headers:
        new_headers.append(new_meta_headers)

    new_headers.append(vcf_reader.column_header)
    file_writer.write("\n".join(new_headers) +"\n")

def _get_execution_metaheaders(pop_values):
    new_meta_headers = []
    for tag, value_list in pop_values.items():
        pop_mean_range, pop_std_range = value_list
        pop_mean_header = "##jacquard.consensus.{0}{1}_RANGE."\
                          "mean_{1}_range={2}"\
                          .format(consensus_helper.JQ_CONSENSUS_TAG, tag,
                                  str(pop_mean_range))
        new_meta_headers.append(pop_mean_header)

        pop_std_header = "##jacquard.consensus.{0}{1}_ZSCORE."\
                         "std_{1}_range={2}"\
                         .format(consensus_helper.JQ_CONSENSUS_TAG, tag,
                                 str(pop_std_range))
        new_meta_headers.append(pop_std_header)

    return "\n".join(new_meta_headers)

def write_to_tmp_file(cons_helper, vcf_reader, tmp_writer):
    vcf_reader.open()
    tmp_writer.open()

    try:
        _write_metaheaders(cons_helper, vcf_reader, tmp_writer)
        logger.info("Adding consensus tags for {}", vcf_reader.input_filepath)
        _add_consensus_tags(cons_helper, vcf_reader, tmp_writer)

    finally:
        vcf_reader.close()
        tmp_writer.close()

def write_to_output_file(cons_helper,
                         execution_context,
                         tmp_reader,
                         file_writer):

    tmp_reader.open()
    file_writer.open()

    try:
        pop_values = cons_helper.get_population_values()

        new_meta_headers = _get_execution_metaheaders(pop_values)
        _write_metaheaders(cons_helper,
                           tmp_reader,
                           file_writer,
                           execution_context,
                           new_meta_headers)

        logger.info("Calculating zscore and writing to {}",
                    file_writer.output_filepath)
        _add_zscore(cons_helper, tmp_reader, file_writer, pop_values)

    finally:
        tmp_reader.close()
        file_writer.close()

def _add_consensus_tags(cons_helper, vcf_reader, file_writer):
    for vcf_record in vcf_reader.vcf_records():
        file_writer.write(cons_helper.add_tags(vcf_record))

def _add_zscore(cons_helper, vcf_reader, file_writer, pop_values):
    for vcf_record in vcf_reader.vcf_records():
        file_writer.write(cons_helper.add_zscore(vcf_record, pop_values))

def add_subparser(subparser):
    # pylint: disable=C0301
    parser = subparser.add_parser("consensus", help="Accepts a Jacquard-merged VCf file and creates a new file, adding consensus fields.")
    parser.add_argument("input", help="Path to Jacquard-merged VCF (or any VCF with Jacquard tags (e.g. JQ_SOM_MT)")
    parser.add_argument("output", help="Path to output VCf")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

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

    cons_helper = consensus_helper.ConsensusHelper()

    vcf_reader = vcf.VcfReader(vcf.FileReader(input_file))
    tmp_output_file = output_file + ".tmp"
    tmp_writer = vcf.FileWriter(tmp_output_file)

    write_to_tmp_file(cons_helper, vcf_reader, tmp_writer)

    tmp_reader = vcf.VcfReader(vcf.FileReader(tmp_output_file))
    file_writer = vcf.FileWriter(output_file)

    write_to_output_file(cons_helper,
                         execution_context,
                         tmp_reader,
                         file_writer)

    os.remove(tmp_output_file)

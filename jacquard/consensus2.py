import numpy as np
import os

from variant_callers import consensus_helper
import vcf as vcf
import utils as utils
import logger as logger

def _write_metaheaders(cons_helper, execution_context,vcf_reader,file_writer):
    new_headers = vcf_reader.metaheaders
    new_headers.extend(execution_context)
    new_headers.extend(cons_helper.get_new_metaheaders())
    new_headers.append(vcf_reader.column_header)
    file_writer.write("\n".join(new_headers) +"\n")
    
def _add_consensus_tags(cons_helper, vcf_reader, file_writer):
    for vcf_record in vcf_reader.vcf_records():
        file_writer.write(cons_helper.add_tags(vcf_record))

def _add_zscore(cons_helper, vcf_reader, file_writer):
    for vcf_record in vcf_reader.vcf_records():
        file_writer.write(cons_helper.add_zscore(vcf_record))

def add_subparser(subparser):
    parser = subparser.add_parser("consensus2", help="Accepts a Jacquard-merged VCf file and creates a new file, adding consensus fields.")
    parser.add_argument("input", help="Path to Jacquard-merged VCF (or any VCF with Jacquard tags (e.g. JQ_SOM_MT)")
    parser.add_argument("output", help="Path to output VCf")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def execute(args, execution_context): 
    input_file = os.path.abspath(args.input)
    output = os.path.abspath(args.output)

    extension = os.path.splitext(os.path.basename(input_file))[1]
    if not os.path.isfile(input_file) or extension != ".vcf":
        raise utils.JQException("Input file [{}] must be a VCF file.", input_file)

    extension = os.path.splitext(os.path.basename(output))[1]
    if extension:
        if extension != ".vcf":
            raise utils.JQException("Output file [{}] must be a VCF file.", output)
        else:
            utils.validate_directories(os.path.dirname(input_file), os.path.dirname(output))
            output_file = output
    else:
        utils.validate_directories(os.path.dirname(input_file), output)
        output_file = os.path.join(output,"consensus.vcf")

    vcf_reader =  vcf.VcfReader(vcf.FileReader(input_file))
    tmp_output_file = output_file + ".tmp"
    tmp_writer = vcf.FileWriter(tmp_output_file)

    cons_helper = consensus_helper.ConsensusHelper()

    vcf_reader.open()
    tmp_writer.open()

    _write_metaheaders(cons_helper, execution_context, vcf_reader, tmp_writer)
    _add_consensus_tags(cons_helper, vcf_reader, tmp_writer)

    vcf_reader.close()
    tmp_writer.close()

#     os.rename(tmp_output_file, output_file)

    tmp_reader = vcf.VcfReader(vcf.FileReader(tmp_output_file))
    file_writer = vcf.FileWriter(output_file)

    tmp_reader.open()
    file_writer.open()

#     _add_zscore(cons_helper, tmp_reader, file_writer)

    tmp_reader.close()
    file_writer.close()

from __future__ import absolute_import
from jacquard.variant_callers import variant_caller_factory
import glob
import os

JQ_OUTPUT_SUFFIX = "jacquardTags"

#TODO: Edit this later to be appropriate for translate.py
def add_subparser(subparser):
    #pylint: disable=line-too-long
    parser = subparser.add_parser("translate", help="Accepts a directory of VCf results and creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to Jacquard-tagged VCFs. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def _mangle_output_filenames(input_file):
    basename, extension = os.path.splitext(os.path.basename(input_file))
    return ".".join([basename, JQ_OUTPUT_SUFFIX, extension.strip(".")])

def _get_output_filenames(input_files):
    output_files = set()
    for input_file in input_files:
        output_files.add(_mangle_output_filenames(input_file))

    return output_files

def _predict_output(args):
    input_dir = os.path.abspath(args.input)
    input_files = sorted(glob.glob(os.path.join(input_dir, "*.vcf")))
    desired_output_files = _get_output_filenames(input_files)

    return desired_output_files

def report_prediction(args):
    return _predict_output(args)

def get_required_input_output_types():
    return ("directory", "directory")

def execute(args, execution_context):
    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)

    in_files = glob.glob(os.path.join(input_dir, "*.vcf"))
#     for in_file in in_files:
#         
#     caller_to_vcf_readers = variant_caller_factory.claim(file_readers)


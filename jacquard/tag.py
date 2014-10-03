from __future__ import print_function
import collections
import glob
import os
import shutil
import sys

from variant_callers import variant_caller_factory
import jacquard_utils
from jacquard_utils import JQException

# pylint: disable=W0142
def _log(msg, *args):
    print(msg.format(*[str(i) for i in args]), file=sys.stderr)


#pylint: disable=C0301
def add_subparser(subparser):
    parser_tag = subparser.add_parser("tag", help="Accepts a directory of VCf results and creates a new directory of VCFs, adding Jacquard-specific FORMAT tags for each VCF record.")
    parser_tag.add_argument("input_dir", help="Path to directory containing VCFs. Other file types ignored")
    parser_tag.add_argument("output_dir", help="Path to Jacquard-tagged VCFs. Will create if doesn't exist and will overwrite files in output directory as necessary")


def tag_files(vcf_readers_to_writers, execution_context):
    total_number_of_files = len(vcf_readers_to_writers)

    for count, in_vcf, vcf_writer in enumerate(vcf_readers_to_writers.items()):
        _log("INFO: Reading [{}] ({}/{})", in_vcf.file_name,
             count, total_number_of_files)
        in_vcf.open()
        vcf_writer.open()

        #TODO cgates: build output header when VcfProvider constructed?
        new_headers = list(in_vcf.metaheaders)
        new_headers.extend(["{}\n".format(header) for header in execution_context])
        new_headers.append("##jacquard.tag.caller={0}\n".format(in_vcf.caller))
        header_block = in_vcf.caller.update_metaheader(''.join(new_headers))
        header_block += in_vcf.column_header

        vcf_writer.write(header_block)

        vcf_writer.write_headers(in_vcf, execution_context)
        for vcf_record in in_vcf.read_record():
            vcf_writer.write(in_vcf.caller.add_tags(vcf_record))

        in_vcf.close()
        vcf_writer.close()


#TODO cgates: add context management to open/close 
class VcfProvider(object):
    def __init__(self, input_dir, file_name, get_caller=None):
        self.input_dir = input_dir
        self.file_name = file_name
        self.name = file_name
        self.file_reader = None
        (self.column_header, self.metaheaders) = self._read_headers()
        self.caller = None
        if get_caller:
            try:
                self.caller = get_caller(self.header())
                _log("DEBUG: VCF [{}] recognized by caller [{}]",
                     self.name, self.caller.name)
            except JQException as ex:
                _log("ERROR: Problem parsing [{}]:{}", self.name, ex)
                raise ex

    def _read_headers(self):
        metaheaders = []
        with open(os.path.join(self.input_dir, self.file_name), "r") as vcf:
            for line in vcf:
                if line.startswith("##"):
                    metaheaders.append(line)
                elif line.startswith("#"):
                    column_header = line
                else:
                    break
        return column_header, metaheaders

    @property
    def header(self):
        return "\n".join(self.metaheaders) + self.column_header + "\n"

    def read_record(self):
        for line in self.file_reader:
            if line.startswith("#"):
                continue
            yield line

    def open(self):
        self.file_reader = open(os.path.join(self.input_dir, self.file_name), 'r')

    def close(self):
        self.file_reader.close()


def _log_caller_info(vcf_providers):
    caller_count = collections.defaultdict(int)
    for vcf in vcf_providers:
        caller_count[vcf.caller.name] += 1
    for caller_name in sorted(caller_count):
        _log("INFO: Recognized [{0}] {1} file(s)",
             caller_count[caller_name], caller_name)


def _build_vcf_providers(in_dir,
                         get_caller=variant_caller_factory.get_caller):
    in_files = sorted(glob.glob(os.path.join(in_dir, "*.vcf")))
    vcf_providers = []
    failures = 0
    for filename in in_files:
        try:
            vcf_providers.append(VcfProvider(in_dir, filename, get_caller))
        except JQException:
            failures += 1
    if failures:
        raise JQException("[{}] VCF files could not be parsed."
                          " Review logs for details, adjust input, "
                          "and try again.", failures)

    _log_caller_info(vcf_providers)

    return vcf_providers


#TODO cgates: add context management to open/close 
class VcfWriter(object):
    
    def __init__(self, output_filepath):
        self.output_filepath = output_filepath
        self.file_writer = None

    def open(self):
        self.file_writer = open(self.output_filepath, "w")

    def write(self, text):
        return self.file_writer.write(text)

    def close(self):
        self.file_writer.close()


def _build_vcf_providers_to_writers(vcf_providers, output_dir):
    vcf_providers_to_writers = {}
    for vcf in vcf_providers:
        basename, extension = os.path.splitext(vcf.file_name)
        new_filename = basename + ".jacquardTags" + extension
        output_filepath = os.path.join(output_dir, new_filename)
        vcf_providers_to_writers[vcf] = VcfWriter(output_filepath)
    return vcf_providers_to_writers

def execute(args, execution_context):
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    jacquard_utils.validate_directories(input_dir, output_dir)

    #TODO cgates: move to jacquard.py
    for line in execution_context:
        _log("DEBUG: {}", line)

    in_vcfs = _build_vcf_providers(input_dir)
    if not in_vcfs:
        _log("ERROR: Specified input directory [{0}] contains no VCF files."
             "Check parameters and try again.", input_dir)
        #TODO cgates: move to jacquard.py
        shutil.rmtree(output_dir)
        exit(1)

    vcfs_to_writers = _build_vcf_providers_to_writers(in_vcfs, output_dir)

    _log("INFO: Processing [{}] VCF file(s) from [{}]", len(in_vcfs), input_dir)
    tag_files(vcfs_to_writers, execution_context)
    _log("INFO: Wrote [{}] VCF file(s) to [{}]",
         len(vcfs_to_writers), output_dir)

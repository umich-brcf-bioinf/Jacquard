#!/usr/bin/env python
##   Copyright 2014 Bioinformatics Core, University of Michigan
##
##   Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##
##       http://www.apache.org/licenses/LICENSE-2.0
##
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

# pylint: disable=missing-docstring,global-statement

from __future__ import absolute_import, print_function

import argparse
import glob
import os
import shutil
import signal
import sys
import traceback

import jacquard.command_validator as command_validator
import jacquard.consensus as consensus
import jacquard.expand as expand
import jacquard.filter_hc_somatic as filter_hc_somatic
import jacquard.logger as logger
import jacquard.merge as merge
import jacquard.merge2 as merge2
import jacquard.normalize as normalize
import jacquard.tag as tag
import jacquard.utils as utils


_SUBCOMMANDS = [normalize,
                tag,
                filter_hc_somatic,
                merge,
                merge2,
                consensus,
                expand]

TMP_DIR_NAME = "jacquard_tmp"
TMP_OUTPUT_PATH = None

def main():
    #pylint: disable=unused-argument
    def handler(signum, frame):
        msg = "WARNING: Jacquard was interrupted before completing."
        try:
            logger.debug(msg)
        finally:
            print(msg, file=sys.stderr)
            exit(1)

    signal.signal(signal.SIGINT, handler)
    signal.signal(signal.SIGTERM, handler)

    dispatch(_SUBCOMMANDS, sys.argv[1:])

def _version_text():
    callers = utils.caller_versions.items()
    caller_versions = [key + " " + value for key, value in callers]
    caller_version_string = "\n\t".join(caller_versions)

    return "Jacquard v{0}\nSupported variant callers:\n\t{1}".\
        format(utils.__version__, caller_version_string)


#TODO: (cgates) This cannot be the simplest thing that could possibly work
def _validate_temp(tmp_output, original_output_dir, force=0):
    extension = os.path.splitext(os.path.basename(tmp_output))[1]
    if extension:
        tmp_dir = os.path.dirname(tmp_output)
    else:
        tmp_dir = tmp_output

    tmp_dir_name = os.path.basename(tmp_dir)
    if os.path.exists(tmp_dir) and not force:
        raise utils.JQException(("A temp directory [{}] already exists in "
                                 "output directory [{}]; move or remove the "
                                 "temp dir or adjust output argument."),
                                tmp_dir_name,
                                original_output_dir)

    if not os.path.exists(tmp_dir):
        try:
            os.mkdir(tmp_dir)
        except:
            message = ("The temp directory [{}] cannot be created in "
                       "output directory [{}]; review permissions "
                       "and output arguments and try again.")
            raise utils.JQException(message,
                                    tmp_dir_name,
                                    original_output_dir)

def _preflight_old(output, desired_output_files, command):
    if os.path.isdir:
        existing_output_paths = sorted(glob.glob(os.path.join(output, "*.vcf")))
    else:
        existing_output_paths = [output]

    existing_output_files = set([os.path.basename(i) for i in existing_output_paths])
    intersection = existing_output_files.intersection(desired_output_files)
    if intersection:
        raise utils.JQException(("ERROR: The command [{}] would "
                                "overwrite existing files {}; review "
                                "command/output dir to avoid overwriting or "
                                "use the flag '--force'. Type 'jacquard -h' "
                                "for more details").format(command,
                                                           list(intersection)))

def _nominate_temp_directory(original_output):
    extension = os.path.splitext(os.path.basename(original_output))[1]

    if extension:
        original_output_fname = os.path.basename(original_output)
        original_output = os.path.dirname(original_output)
        tmp_output = os.path.join(original_output,
                                  TMP_DIR_NAME,
                                  original_output_fname)
    else:
        tmp_output = os.path.join(original_output, TMP_DIR_NAME)

    return tmp_output

def _create_temp_directory(original_output_dir, force=0):
    extension = os.path.splitext(os.path.basename(original_output_dir))[1]

    if extension:
        original_output_fname = os.path.basename(original_output_dir)
        original_output_dir = os.path.dirname(original_output_dir)
        tmp_output = os.path.join(original_output_dir,
                                  TMP_DIR_NAME,
                                  original_output_fname)
    else:
        tmp_output = os.path.join(original_output_dir, TMP_DIR_NAME)

    try:
        os.mkdir(original_output_dir)
    except OSError:
        pass

    if len(os.listdir(original_output_dir)) != 0:
        if not force:
            raise utils.JQException("Specified output {} is not empty. "
                                    "Please specify an empty directory or "
                                    "use the flag '--force'. Type "
                                    "'jacquard -h' for more details.",
                                    original_output_dir)

    _validate_temp(tmp_output, original_output_dir, force)

    return tmp_output

def _move_tmp_contents_to_original(tmp_path, original_output):
    if os.path.isfile(tmp_path):
        tmp_path = os.path.dirname(tmp_path)

    tmp_contents = os.listdir(tmp_path)
    for fname in tmp_contents:
        full_fname = os.path.join(tmp_path, fname)

        if os.path.isfile(full_fname):
            shutil.move(full_fname, original_output)
        else:
            if not os.path.isdir(original_output):
                output_file = os.path.dirname(original_output)
                shutil.move(full_fname, output_file)

    os.rmdir(tmp_path)

def dispatch(modules, arguments):
    # pylint: disable=line-too-long
    parser = argparse.ArgumentParser(
        usage="jacquard",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''type 'Jacquard -h <subcommand>' for help on a specific subcommand''',
        epilog="authors: Jessica Bene, Ashwini Bhasi, Chris Gates, Kevin Meng, Peter Ulintz; October 2014")

    parser.add_argument("-V",
                        "--version",
                        action='version',
                        version=_version_text())

    subparsers = parser.add_subparsers(title="subcommands",
                                       dest="subparser_name")

    try:
        module_dispatch = {}
        for module in modules:
            module.add_subparser(subparsers)
            short_name = module.__name__.split('.')[-1]
            module_dispatch[short_name] = module

        cwd = os.path.dirname(os.getcwd())
        execution_context = [\
            "##jacquard.version={0}".format(utils.__version__),
            "##jacquard.command={0}".format(" ".join(arguments)),
            "##jacquard.cwd={0}".format(cwd)]

        args = parser.parse_args(arguments)

        logger.initialize_logger(args.subparser_name)
        logger.info("Jacquard begins (v{})", utils.__version__)
        logger.info("Saving log to [{}]", logger.log_filename)
        logger.debug("cwd|{}", os.getcwd())
        logger.debug("command|{}", " ".join(arguments))

        original_output = args.output

        command_validator.preflight(args, module_dispatch[args.subparser_name])

        TMP_OUTPUT_PATH = _nominate_temp_directory(original_output)
        args.output = TMP_OUTPUT_PATH
        logger.debug("Writing output to tmp directory [{}]", TMP_OUTPUT_PATH)

        module_dispatch[args.subparser_name].execute(args, execution_context)

        logger.debug("Moving files from tmp directory {} to output directory",
                     TMP_OUTPUT_PATH,
                     original_output)
        _move_tmp_contents_to_original(TMP_OUTPUT_PATH, original_output)
        logger.debug("Removed tmp directory {}", TMP_OUTPUT_PATH)

        logger.debug("Output saved to [{}]", os.path.abspath(original_output))
        logger.info("Output saved to [{}]", original_output)
        if logger.SHOW_WARNING:
            logger.info("Done. (See warnings above)")
        else:
            logger.info("Done")

    #pylint: disable=broad-except
    except Exception as exception:
        logger.error(str(exception))
        logger.error("Jacquard encountered an unanticipated problem. Please review log file and contact your sysadmin or Jacquard support for assistance.")
        logger.debug(traceback.format_exc())
        sys.exit(1)
    finally:
        _cleanup()


def _cleanup():
    if TMP_OUTPUT_PATH and os.path.exists(TMP_OUTPUT_PATH):
        logger.debug("Cleaning up tmp directory")
        shutil.rmtree(TMP_OUTPUT_PATH)

if __name__ == '__main__':
    main()

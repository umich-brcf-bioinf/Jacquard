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

# pylint: disable=C0111
# pylint: disable-msg=W0403

from __future__ import print_function
import argparse
import os
import signal
import shutil
import sys
import traceback

import tag as tag
import normalize as normalize
import filter_hc_somatic as filter_hc_somatic
import merge as merge
import consensus as consensus
import expand
import utils as utils
import logger as logger


_SUBCOMMANDS = [normalize,
                tag,
                filter_hc_somatic,
                merge,
                consensus,
                expand]

TMP_DIR_NAME = "jacquard.tmp"
TMP_OUTPUT_PATH = None

def main():
    
    #pylint: disable=W0613
    def handler(signum, frame):
        print("WARNING: Jacquard was interrupted before completing.",
              file=sys.stderr)
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


#TODO (cgates): This cannot be the simplest thing that could possibly work
def _validate_temp(tmp_output, original_output_dir, force=0):
    if os.path.isfile(original_output_dir):
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

def _create_temp_directory(original_output_dir, force=0):
    if os.path.isfile(original_output_dir):
        original_output_fname = os.path.basename(original_output_dir)
        original_output_dir = os.path.dirname(original_output_dir)
        tmp_output = os.path.join(original_output_dir,
                                  TMP_DIR_NAME,
                                  original_output_fname)
    else:
        tmp_output = os.path.join(original_output_dir, TMP_DIR_NAME)

    try:
        os.mkdir(original_output_dir)
    except:
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

def _move_tmp_contents_to_original(tmp_dir, original_output):
    if os.path.isfile(tmp_dir):
        tmp_dir = os.path.dirname(tmp_dir)

    tmp_contents = os.listdir(tmp_dir)
    for fname in tmp_contents:
        full_fname = os.path.join(tmp_dir, fname)

        if os.path.isfile(full_fname):
            shutil.move(full_fname, original_output)
        else:
            if not os.path.isdir(original_output):
                output_dir = os.path.dirname(original_output)
                shutil.move(full_fname, output_dir)

    os.rmdir(tmp_dir)

# pylint: disable=C0301
def dispatch(modules, arguments):
    
    parser = argparse.ArgumentParser(
        usage="jacquard",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''type 'Jacquard -h <subcommand>' for help on a specific subcommand''',
        epilog="authors: Jessica Bene, Ashwini Bhasi, Chris Gates, Kevin Meng, Peter Ulintz; October 2014")

    parser.add_argument(\
                        "-V",
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

        original_output_dir = args.output
#
#        global TMP_OUTPUT_PATH
#        TMP_OUTPUT_PATH = _create_temp_directory(original_output_dir, args.force)
#        args.output = TMP_OUTPUT_PATH
#        logger.debug("Writing output to tmp directory [{}]", TMP_OUTPUT_PATH)

        module_dispatch[args.subparser_name].execute(args, execution_context)

#        logger.debug("Moving files from tmp directory {} to output directory", TMP_OUTPUT_PATH, original_output_dir)
#        _move_tmp_contents_to_original(TMP_OUTPUT_PATH, original_output_dir)
#        logger.debug("Removed tmp directory {}", TMP_OUTPUT_PATH)

        logger.info("Output saved to [{}]", original_output_dir)
        logger.info("Done")

    # pylint: disable=W0703
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

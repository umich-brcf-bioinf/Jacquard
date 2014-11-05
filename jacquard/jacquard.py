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

from __future__ import print_function
import argparse
import os
import signal
import sys
import traceback

import tag as tag
import normalize as normalize
import filter_hc_somatic as filter_hc_somatic
import merge as merge
import consensus as consensus
import expand_old
import expand
import utils as utils

import logger


_SUBCOMMANDS=[normalize,
              tag,
              filter_hc_somatic,
              merge,
              consensus,
              expand_old,
              expand]

def main():
    #pylint: disable=W0613
    def handler(signum, frame):
        print("WARNING: Jacquard was interrupted before completing.",
              file=sys.stderr)
        exit(1)

    signal.signal(signal.SIGINT, handler)
    signal.signal(signal.SIGTERM, handler)

    dispatch(_SUBCOMMANDS, sys.argv[1:])

def version_text():
    callers = utils.caller_versions.items()
    caller_versions = [key + " " + value for key, value in callers]
    caller_version_string = "\n\t".join(caller_versions)
    
    return "Jacquard v{0}\nSupported variant callers:\n\t{1}".\
        format(utils.__version__, caller_version_string)

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
                        version=version_text())

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


        module_dispatch[args.subparser_name].execute(args, execution_context)

    # pylint: disable=W0703
    except Exception as exception:
        logger.error(str(exception))
        logger.error("Jacquard encountered an unanticipated problem. Please review log file and contact your sysadmin or Jacquard support for assistance.")
        logger.debug(traceback.format_exc())
        sys.exit(1)

    logger.info("Done")

if __name__ == '__main__':
    main()

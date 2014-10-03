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

import pivot_variants
import rollup_genes
import tag
import normalize
import filter_hc_somatic
import merge
import consensus
import expand
import style
import jacquard_utils

_SUBCOMMANDS=[normalize,
              tag,
              filter_hc_somatic,
              merge,
              consensus,
              expand,
              style,
              rollup_genes,
              pivot_variants]

def main():
    def handler(signum, frame):
        print("WARNING: Jacquard was interrupted before completing.", file=sys.stderr)
        exit(1)

    signal.signal(signal.SIGINT, handler)
    signal.signal(signal.SIGTERM, handler)

    dispatch(_SUBCOMMANDS, sys.argv[1:])

def version_text():
    callers = jacquard_utils.caller_versions.items()
    caller_versions = [key + " " + value for key, value in callers]
    caller_version_string = "\n\t".join(caller_versions)
    return "Jacquard v{0}\nSupported variant callers:\n\t{1}".\
        format(jacquard_utils.__version__, caller_version_string)

# pylint: disable=C0301
def dispatch(modules, arguments):
    parser = argparse.ArgumentParser(
        usage="jacquard",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''type 'Jacquard -h <subcommand>' for help on a specific subcommand''',
        epilog="authors: Jessica Bene, Chris Gates 07/2014")
    parser.add_argument(\
                        "-v",
                        "--version",
                        action='version',
                        version=version_text())
    subparsers = parser.add_subparsers(title="subcommands",
                                       dest="subparser_name")

    module_dispatch = {}
    for module in modules:
        module.add_subparser(subparsers)
        short_name = module.__name__.split('.')[-1]
        module_dispatch[short_name] = module

    execution_context = [\
        "##jacquard.version={0}".format(jacquard_utils.__version__),
        "##jacquard.command={0}".format(" ".join(arguments)),
        "##jacquard.cwd={0}".format(os.path.dirname(os.getcwd()))]

    args = parser.parse_args(arguments)
    module_dispatch[args.subparser_name].execute(args, execution_context)

if __name__ == '__main__':
    main()

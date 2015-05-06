#!/usr/bin/env python
"""Entry-point for suite of sub-commands.

The only executable module in the project; this module
 * validates command line args
 * manages use of temp directories (to keep output clean and atomic)
 * dipatches to sub-commands as appropriate
 * attempts to deal with usage and run-time errors

Jacquard first writes results to temp dir and only copies results on successful
    completion.

Then architecture of Jacquard modules can be divided into:
 * commands : These transform files or directories (e.g. translate.py) and are
     indirectly executable through the jacquard module. Each command must
     implement an execute method that does the heavy lifting along with some
     simpler methods that expedite command validation
 * callers : These transform VcfRecords (e.g. mutect.py). They typically have
     a collection of tag classes; where each tag holds the metaheader and
     code to transform a single VcfRecord. Note that a caller could
     manipulate any aspect of a VcfRecord, but (by strong convention) typically
     only adds information, for example add a sample-format tag, add an info
     field, or add a filter field.
 * helpers : Common functionality (e.g. command_validator, vcf, logger, etc.)
"""
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
from __future__ import print_function, absolute_import, division

import argparse
from datetime import datetime
import distutils.dir_util
import os
import re
import shutil
import signal
import sys
import traceback

from jacquard import __version__
import jacquard.utils.command_validator as command_validator
import jacquard.expand as expand
import jacquard.utils.logger as logger
import jacquard.merge as merge
import jacquard.summarize as summarize
import jacquard.translate as translate
import jacquard.utils.utils as utils
from jacquard.variant_caller_transforms import variant_caller_factory


_SUBCOMMANDS = [translate,
                merge,
                summarize,
                expand]


class _JacquardArgumentParser(argparse.ArgumentParser):
    """Argument parser that raises UsageError instead of exiting."""
    #pylint: disable=too-few-public-methods

    @staticmethod
    def _remessage_invalid_subparser(message):
        pattern = r"argument subparser_name: invalid choice: '(.*?)' \((.*)\)"
        match = re.match(pattern, message)
        if match:
            message = ("'{}' is not a Jacquard command;"
                       " {}.").format(match.group(1), match.group(2))
        return message

    def error(self, message):
        '''Suppress default exit behavior'''
        message = self._remessage_invalid_subparser(message)
        raise utils.UsageError(message)


def _cleanup(temp_working_dir):
    if temp_working_dir and os.path.exists(temp_working_dir):
        logger.debug("Cleaning up tmp directory")
        shutil.rmtree(temp_working_dir)


def _move_tmp_contents_to_original(args):
    temp_dir = args.temp_working_dir
    dest_dir = os.path.dirname(args.temp_working_dir)

    logger.debug("Moving files from tmp directory [{}] to output "
                 "directory [{}]",
                 temp_dir,
                 dest_dir)
    distutils.dir_util.copy_tree(temp_dir, dest_dir)
    logger.debug("Output saved to [{}]", os.path.abspath(args.original_output))
    logger.info("Output saved to [{}]", args.original_output)
    _cleanup(temp_dir)


def _parse_command_line_args(modules, arguments):
    parser = _JacquardArgumentParser(
        usage="jacquard",
#         formatter_class=argparse.RawDescriptionHelpFormatter,
        formatter_class=argparse.RawTextHelpFormatter,
        # pylint: disable=line-too-long
        description='''type 'jacquard <command> -h' for help on a specific command''',
        epilog="See https://github.com/umich-brcf-bioinf/Jacquard for more info.")

    parser.add_argument("-V",
                        "--version",
                        action='version',
                        version=_version_text())
    subparsers = parser.add_subparsers(title="subcommands",
                                       dest="subparser_name")

    module_dispatch = {}
    for module in modules:
        module.add_subparser(subparsers)
        short_name = module.__name__.split('.')[-1]
        module_dispatch[short_name] = module

    args = parser.parse_args(arguments)
    return module_dispatch[args.subparser_name], args


def _set_interrupt_handler(target=signal.signal):
    def _handler(dummy_signum, dummy_frame):
        msg = "WARNING: Jacquard was interrupted before completing."
        try:
            logger.debug(msg)
        finally:
            print(msg, file=sys.stderr)
            exit(1)
    target(signal.SIGINT, _handler)
    target(signal.SIGTERM, _handler)


def _version_text():
    callers = variant_caller_factory.SUPPORTED_CALLER_VERSIONS.items()
    caller_versions = [key + " " + value for key, value in callers]
    caller_version_string = "\n\t".join(caller_versions)
    return "Jacquard v{0}\nSupported variant callers:\n\t{1}".\
        format(__version__, caller_version_string)

def _get_execution_context(command):
    cwd = os.path.dirname(os.getcwd())
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    return ['##jacquard=<Timestamp="{}",Command="{}",Cwd="{}">'.format(now,
                                                                       command,
                                                                       cwd)]

def _dispatch(modules, arguments):
    try:
        command, args = _parse_command_line_args(modules, arguments)
        execution_context = _get_execution_context(command)

        logger.initialize_logger(args)
        logger.debug("Jacquard run begins")
        logger.debug("cwd|{}", os.getcwd())
        logger.debug("command|{}", " ".join(arguments))

        command_validator.preflight(command, args)

        logger.info("Jacquard begins (v{})", __version__)
        logger.info("Saving log to [{}]",
                    os.path.basename(logger.log_filename))
        logger.debug("Writing output to tmp directory [{}]",
                     args.temp_working_dir)

        command.execute(args, execution_context)

        _move_tmp_contents_to_original(args)

        if logger.WARNING_OCCURRED:
            logger.info("Done. (See warnings above)")
        else:
            logger.info("Done")

    except utils.UsageError as usage_error:
        message = "Jacquard usage problem: {}".format(str(usage_error))
        logger.debug(message)
        print(message, file=sys.stderr)
        try:
            print("See 'jacquard {} --help'.".format(args.subparser_name),
                  file=sys.stderr)
        except NameError: #could not determine the command
            print("See 'jacquard --help'.", file=sys.stderr)
        sys.exit(1)
    except Exception as exception: #pylint: disable=broad-except
        logger.error(str(exception))
        logger.error("Jacquard encountered an unanticipated problem. "
                     "Please review log file and contact your sysadmin "
                     "or Jacquard support for assistance.")
        logger.debug(traceback.format_exc())
        sys.exit(1)
    finally:
        try:
            _cleanup(args.temp_working_dir)
        except Exception: #pylint: disable=broad-except
            pass #we tried


def main():
    _set_interrupt_handler()
    _dispatch(_SUBCOMMANDS, sys.argv[1:])


if __name__ == '__main__':
    main()

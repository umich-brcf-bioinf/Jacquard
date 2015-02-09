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


from __future__ import absolute_import, print_function

import argparse
import distutils.dir_util
import os
import re
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
import jacquard.normalize as normalize
import jacquard.tag as tag
import jacquard.utils as utils


_SUBCOMMANDS = [normalize,
                tag,
                filter_hc_somatic,
                merge,
                consensus,
                expand]


class _JacquardArgumentParser(argparse.ArgumentParser):
    '''Suppress default exit behavior'''

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
        formatter_class=argparse.RawDescriptionHelpFormatter,
        # pylint: disable=line-too-long
        description='''type 'Jacquard -h <command>' for help on a specific command''',
        epilog="authors: Jessica Bene, Ashwini Bhasi, Chris Gates, Kevin Meng, Peter Ulintz; October 2014")

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
    def _handler(signum, frame): #pylint: disable=unused-argument
        msg = "WARNING: Jacquard was interrupted before completing."
        try:
            logger.debug(msg)
        finally:
            print(msg, file=sys.stderr)
            exit(1)
    target(signal.SIGINT, _handler)
    target(signal.SIGTERM, _handler)


def _version_text():
    callers = utils.caller_versions.items()
    caller_versions = [key + " " + value for key, value in callers]
    caller_version_string = "\n\t".join(caller_versions)
    return "Jacquard v{0}\nSupported variant callers:\n\t{1}".\
        format(utils.__version__, caller_version_string)


def dispatch(modules, arguments):
    try:
        command, args = _parse_command_line_args(modules, arguments)

        cwd = os.path.dirname(os.getcwd())
        execution_context = [\
            "##jacquard.version={0}".format(utils.__version__),
            "##jacquard.command={0}".format(" ".join(arguments)),
            "##jacquard.cwd={0}".format(cwd)]

        logger.initialize_logger(args.subparser_name)
        logger.debug("Jacquard run begins")
        logger.debug("cwd|{}", os.getcwd())
        logger.debug("command|{}", " ".join(arguments))

        command_validator.preflight(args, command)

        logger.info("Jacquard begins (v{})", utils.__version__)
        logger.info("Saving log to [{}]", logger.log_filename)
        logger.debug("Writing output to tmp directory [{}]",
                     args.temp_working_dir)

        command.execute(args, execution_context)

        _move_tmp_contents_to_original(args)

        if logger.SHOW_WARNING:
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
    dispatch(_SUBCOMMANDS, sys.argv[1:])


if __name__ == '__main__':
    main()

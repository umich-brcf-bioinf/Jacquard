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
import os
import re
import shutil
import signal
import sys
import time
import traceback


_SUBCOMMANDS = [normalize,
                tag,
                filter_hc_somatic,
                merge,
                merge2,
                consensus,
                expand]

_TEMP_WORKING_DIR_FORMAT = "jacquard.{}.{}.tmp"


class _JacquardArgumentParser(argparse.ArgumentParser):
    '''Suppress default exit behavior'''
    def error(self, message):
        message = self._remessage_invalid_subparser(message)
        raise utils.UsageError(message)

    @staticmethod
    def _remessage_invalid_subparser(message):
        pattern = r"argument subparser_name: invalid choice: '(.*?)' \((.*)\)"
        match = re.match(pattern, message)
        if match:
            message = ("'{}' is not a Jacquard command;"
                       " {}.").format(match.group(1), match.group(2))
        return message


def main():
    _set_interrupt_handler()
    dispatch(_SUBCOMMANDS, sys.argv[1:])


def _set_interrupt_handler(target=signal.signal):
    def handler(signum, frame): #pylint: disable=unused-argument
        msg = "WARNING: Jacquard was interrupted before completing."
        try:
            logger.debug(msg)
        finally:
            print(msg, file=sys.stderr)
            exit(1)

    target(signal.SIGINT, handler)
    target(signal.SIGTERM, handler)


def _version_text():
    callers = utils.caller_versions.items()
    caller_versions = [key + " " + value for key, value in callers]
    caller_version_string = "\n\t".join(caller_versions)

    return "Jacquard v{0}\nSupported variant callers:\n\t{1}".\
        format(utils.__version__, caller_version_string)

def _get_temp_working_dir(original_output):
    base_dir = os.path.dirname(os.path.abspath(original_output))
    pid = os.getpid()
    microseconds_since_epoch = int(time.time() * 1000 * 1000)
    dir_name = _TEMP_WORKING_DIR_FORMAT.format(str(pid),
                                               str(microseconds_since_epoch))
    temp_working_dir = os.path.join(base_dir, dir_name)
    new_output = os.path.join(temp_working_dir,
                              os.path.basename(original_output))

    return temp_working_dir, new_output

def _move_tmp_contents_to_original(args):
    temp_dir = args.temp_working_dir
    dest_dir = os.path.dirname(args.temp_working_dir)

    logger.debug("Moving files from tmp directory {} to output directory",
                 temp_dir,
                 dest_dir)
    distutils.dir_util.copy_tree(temp_dir, dest_dir)
    _cleanup(temp_dir)

def dispatch(modules, arguments):
    try:
        # pylint: disable=line-too-long
        parser = _JacquardArgumentParser(
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

        module_dispatch = {}
        for module in modules:
            module.add_subparser(subparsers)
            short_name = module.__name__.split('.')[-1]
            module_dispatch[short_name] = module

        args = parser.parse_args(arguments)

        cwd = os.path.dirname(os.getcwd())
        execution_context = [\
            "##jacquard.version={0}".format(utils.__version__),
            "##jacquard.command={0}".format(" ".join(arguments)),
            "##jacquard.cwd={0}".format(cwd)]

        logger.initialize_logger(args.subparser_name)
        logger.debug("Jacquard run begins")
        logger.debug("cwd|{}", os.getcwd())
        logger.debug("command|{}", " ".join(arguments))

        args.original_output = args.output
        (args.temp_working_dir,\
            args.output) = _get_temp_working_dir(args.original_output)

        command_validator.preflight(args, module_dispatch[args.subparser_name])

        logger.info("Jacquard begins (v{})", utils.__version__)
        logger.info("Saving log to [{}]", logger.log_filename)
        logger.debug("Writing output to tmp directory [{}]",
                     args.temp_working_dir)

        module_dispatch[args.subparser_name].execute(args, execution_context)

        _move_tmp_contents_to_original(args)

        logger.debug("Output saved to [{}]",
                     os.path.abspath(args.original_output))
        logger.info("Output saved to [{}]",
                    args.original_output)
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
        except NameError:
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


def _cleanup(temp_working_dir):
    if temp_working_dir and os.path.exists(temp_working_dir):
        logger.debug("Cleaning up tmp directory")
        shutil.rmtree(temp_working_dir)

if __name__ == '__main__':
    main()

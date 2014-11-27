#!/usr/bin/env python2.7
from __future__ import print_function
from datetime import datetime
import logging
import getpass
import os
import socket
import sys

SHOW_WARNING = False

_DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
_FILE_LOG_FORMAT = ('%(asctime)s|%(levelname)s|%(start_time)s|%(host)s|%(user)s'
                    '|%(tool)s|%(message)s')
_CONSOLE_LOG_FORMAT = '%(asctime)s|%(levelname)s|%(tool)s|%(message)s'
# pylint: disable=C0103
logging_dict = {}
_verbose = False
log_filename = None

#pylint: disable=W0603
def initialize_logger(tool, verbose=False):
    global log_filename
    log_filename = os.path.join(os.getcwd(), "jacquard.log")
    logging.basicConfig(format=_FILE_LOG_FORMAT,
                        level="DEBUG",
                        datefmt=_DATE_FORMAT,
                        filename=log_filename)

    global _verbose
    _verbose = verbose

    start_time = datetime.now().strftime(_DATE_FORMAT)
    global logging_dict
    logging_dict = {'user': getpass.getuser(),
                    'host': socket.gethostname(),
                    'start_time': start_time,
                    'tool': tool}

def error(message, *args):
    _print("ERROR", message, args)
    logging.error(_format(message, args), extra=logging_dict)

def warning(message, *args):
    _print("WARNING", message, args)
    logging.warning(_format(message, args), extra=logging_dict)
    global SHOW_WARNING
    SHOW_WARNING = True

def info(message, *args):
    _print("INFO", message, args)
    logging.info(_format(message, args), extra=logging_dict)

def debug(message, *args):
    if _verbose:
        _print("DEBUG", message, args)
    logging.debug(_format(message, args), extra=logging_dict)

def _print(level, message, args):
    now = datetime.now().strftime(_DATE_FORMAT)
    print(_CONSOLE_LOG_FORMAT % {'asctime': now,
                                 'levelname':level,
                                 'tool':logging_dict['tool'],
                                 'message': _format(message, args)},
          file=sys.stderr)

# pylint: disable=W0703
def _format(message, args):
    try:
        log_message = message.format(*[str(i) for i in args])
    except Exception as err:
        log_message = ("Malformed log message ({}: {})"
                       "|{}|{}").format(type(err).__name__,
                                        err,
                                        message,
                                        [str(i) for i in args])
    return log_message

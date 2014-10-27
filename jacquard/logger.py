#!/usr/bin/env python2.7
from __future__ import print_function
from datetime import datetime
import logging
import getpass
import os
import shutil
import socket
import sys

_FILE_LOG_FORMAT = '%(asctime)s|%(levelname)s|%(time)s|%(host)s|%(user)s|%(tool)s|%(message)s'  ### File logs are not done yet
_CONSOLE_LOG_FORMAT = '%(asctime)s|%(levelname)s|%(tool)s|%(message)s' ### Console prints are tested.
logging_dict = {}

def initialize_logger(tool):
    log_dir=os.path.join(os.path.dirname(os.getcwd()), "logs")
    if not os.path.isdir(log_dir):
        os.mkdir(log_dir)
    time = datetime.now()
    logging.basicConfig(format=_FILE_LOG_FORMAT, level="DEBUG", datefmt='%Y/%m/%d %I:%M:%S %p', filename=os.path.join(log_dir, "jacquard.log"))
    
    global logging_dict
    logging_dict = {'user': getpass.getuser(), 'host': socket.gethostname(), 'time': time, 'tool': tool}
    
# def error(message, logging_dict = {}, tool = ""):
def error(message, *args):
    _printer("ERROR", message, *args)
    logging.error(message.format(*[str(i) for i in args]), extra=logging_dict)

# def warning(message, logging_dict = {}, tool = ""):
def warning(message, *args):
    _printer("WARNING", message, *args)
    logging.warning(message.format(*[str(i) for i in args]), extra=logging_dict)
    
def info(message, *args):
    _printer("INFO", message, *args)
    logging.info(message.format(*[str(i) for i in args]), extra=logging_dict)

def debug(message, *args):# def debug(message, logging_dict = {}, tool = ""):
#     _printer("DEBUG", message, *args)
    logging.debug(message.format(*[str(i) for i in args]), extra=logging_dict)

def _printer(level, message, *args):
    print (_CONSOLE_LOG_FORMAT % {'asctime':datetime.now().strftime('%Y/%m/%d %I:%M:%S %p'),
                                  'levelname':level, 
                                  'tool':logging_dict['tool'], 
                                  'message':message.format(*[str(i) for i in args])}, file=sys.stderr)
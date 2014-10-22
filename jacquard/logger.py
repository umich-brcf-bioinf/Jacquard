#!/usr/bin/env python2.7
from datetime import datetime
import logging
import getpass
import os
import shutil
import socket

def initialize_logger(output_dir, tool):
    log_dir =os.path.join(output_dir + "/logs/")
    if not os.path.isdir(log_dir):
        os.mkdir(log_dir)
    time = datetime.now()
    logging.basicConfig(format='%(asctime)s|%(levelname)s|%(time)s|%(host)s|%(user)s|%(tool)s|%(message)s', level="DEBUG", datefmt='%Y/%m/%d %I:%M:%S %p', filename= output_dir + "/logs/jacquard.log")
   
    logging_dict = {'user': getpass.getuser(), 'host': socket.gethostname(), 'time': time,'tool': tool}

    return logging_dict
    
def log_error(message, logging_dict = {}, tool = ""):
    logging.error(message, extra=logging_dict)
    
def log_warning(message, logging_dict = {}, tool = ""):
    print "Warning: " + message
    logging.warning(message, extra=logging_dict)
    
def log_info(message, logging_dict = {}, tool = ""):
    print "{}|{}|{}".format(datetime.now(), tool, message)
    logging.info(message, extra=logging_dict)

def log_debug(message, logging_dict = {}, tool = ""):
    logging.debug(message, extra=logging_dict)

def log_validation_messages(message, warning_message, logging_dict):
    if warning_message != "":
        log_warning(warning_message, logging_dict)
    if message != "":
        log_error(message, logging_dict)

    
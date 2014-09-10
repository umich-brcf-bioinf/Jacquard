#!/usr/bin/env python2.7
import logging
import getpass
import os
import shutil
import socket

def initialize_logger(output_dir):
    log_dir =os.path.join(output_dir + "/logs/")
    if not os.path.isdir(log_dir):
        os.mkdir(log_dir)
        
    logging.basicConfig(format='%(asctime)s|%(user)s|%(host)s|%(levelname)s|%(message)s', level="DEBUG", datefmt='%Y/%m/%d %I:%M:%S %p', filename= output_dir + "/logs/jacquard.log")
   
    logging_dict = {'user': getpass.getuser(), 'host': socket.gethostname()}

    return logging_dict
    
def log_error(message, logging_dict = {}):
    print "Error: " + message
    logging.error(message, extra=logging_dict)
    exit(1)
    
def log_warning(message, logging_dict = {}):
    print "Warning: " + message
    logging.warning(message, extra=logging_dict)
    
def log_success(message, logging_dict = {}):
    print message
    logging.info(message, extra=logging_dict)

def log_validation_messages(message, warning_message, logging_dict):
    if warning_message != "":
        log_warning(warning_message, logging_dict)
    if message != "":
        log_error(message, logging_dict)

    
#!/usr/bin/python2.7
import argparse

execute_called = False

def add_subparser(subparser):
    parser_mock_module = subparser.add_parser("test.mock_module", help="foo")
    
def execute(args, execution_context):
    global execute_called
    execute_called = True

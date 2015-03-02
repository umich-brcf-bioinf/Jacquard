#pylint: disable=invalid-name, global-statement
from __future__ import print_function, absolute_import
from collections import defaultdict

messages = defaultdict(list)

def reset():
    global messages
    messages = defaultdict(list)

def error(message, *args):
    messages["ERROR"].append(_format(message, args))

def warning(message, *args):
    messages["WARNING"].append(_format(message, args))

def info(message, *args):
    messages["INFO"].append(_format(message, args))

def debug(message, *args):
    messages["DEBUG"].append(_format(message, args))

def _format(message, args):
    return message.format(*[str(i) for i in args])


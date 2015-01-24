#pylint: disable=invalid-name,too-few-public-methods
from __future__ import absolute_import, print_function
from os import listdir
import collections
import jacquard.logger as logger
import os
import copy
import re
#TODO cgates: This should be in jacquard.__init__?
__version__ = 0.21

#TODO: cgates: These should be in the callers or the caller factory, but not here.
global caller_versions
caller_versions = {"VarScan":"v2.3", "MuTect": "v1.1.4", "Strelka": "v2.0.15"}

global jq_somatic_tag
global jq_af_tag
global jq_dp_tag
jq_somatic_tag = "HC_SOM"
jq_af_tag = "AF"
jq_dp_tag = "DP"

class JQException(Exception):
    """Base class for exceptions in this module."""
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(JQException, self).__init__(error_msg)

#TODO: (cgates): Suspect this should raise exception instead of logging and exiting?
def validate_directories(input_dir=None, output_dir=None):
    if input_dir:
        if not os.path.isdir(input_dir):
            logger.error("Specified input directory [{}] does not exist.",
                         input_dir)
            exit(1)
        try:
            listdir(input_dir)
        except OSError:
            logger.error("Specified input directory [{}] cannot be read. "+
                         "Check permissions and try again.", input_dir)
            exit(1)
    if output_dir:
        if not os.path.isdir(output_dir):
            try:
                os.makedirs(output_dir)
            except OSError:
                logger.error("Output directory [{}] could not be created. "+
                             "Check parameters and try again", output_dir)
                exit(1)

class OrderedSet(collections.MutableSet):
    def __init__(self, iterable=None):
        self.end = end = []
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:
            key, prev, nxt = self.map.pop(key)
            prev[2] = nxt
            nxt[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    #pylint: disable=arguments-differ
    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)

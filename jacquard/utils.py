#pylint: disable=invalid-name,too-few-public-methods, global-at-module-level
from __future__ import absolute_import, print_function
import collections

#TODO: cgates: These should be in the callers/caller factory, but not here.
global caller_versions
caller_versions = {"VarScan":"v2.3", "MuTect": "v1.1.4", "Strelka": "v2.0.15"}

global jq_somatic_tag
global jq_af_tag
global jq_dp_tag
jq_somatic_tag = "HC_SOM"
jq_af_tag = "AF"
jq_dp_tag = "DP"

def round_two_digits(val):
    if len(val.split(".")[1]) > 2:
        return str(round(100 * float(val))/100)
    return val
#     return "{0:.2f}".format(float(val))


#TODO: (jebene/kmeng) - do we want to truncate decimals if it's X.0?
#     if len(val.split(".")[1]) <= 2:
#         if val.split(".")[1] == '0':
#             return val.split(".")[0]
#         return val
#
#     else:
#         return str(round(100 * float(val))/100)

class JQException(Exception):
    """Base class for exceptions in this module."""
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(JQException, self).__init__(error_msg)


class UsageError(JQException):
    """Raised for malformed command or invalid arguments"""
    def __init__(self, msg, *args):
        super(UsageError, self).__init__(msg, *args)


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

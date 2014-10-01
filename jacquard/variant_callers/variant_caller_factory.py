from jacquard.variant_callers.varscan import Varscan
from jacquard.variant_callers.strelka import Strelka
from jacquard.variant_callers.mutect import Mutect

from jacquard.jacquard_utils import JQException

class JQCallerNotRecognizedException(JQException):
    pass

class Factory():
    def __init__(self):
        self._callers = [Varscan(), Strelka(), Mutect()]
        
    def get_caller(self, header):
        for caller in self._callers:
            if caller.validate_input_file(header):
                return caller
        raise JQCallerNotRecognizedException()
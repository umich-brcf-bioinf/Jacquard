from jacquard.variant_callers.varscan import Varscan
from jacquard.variant_callers.strelka import Strelka
from jacquard.variant_callers.mutect import Mutect
from jacquard.variant_callers.unknown import Unknown

from jacquard.jacquard_utils import JQException

class Factory():
    def __init__(self):
        self._callers = [Varscan(), Strelka(), Mutect(), Unknown()]
        
    def get_caller(self, header):
        for caller in self._callers:
            if caller.validate_input_file(header):
                if caller.name == "Unknown":
                    raise JQException
                else:
                    return caller
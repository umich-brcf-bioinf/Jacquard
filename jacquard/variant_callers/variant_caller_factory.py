from jacquard.variant_callers.varscan import Varscan
from jacquard.variant_callers.strelka import Strelka
from jacquard.variant_callers.mutect import Mutect

from jacquard.jacquard_utils import JQException

#TODO cgates: suspect this does not need to be a class
class Factory(object):
    def __init__(self):
        self._callers = [Varscan(), Strelka(), Mutect()]

    def get_caller(self, vcf):
        for caller in self._callers:
            if caller.validate_input_file(vcf.header):
                return caller
        raise JQException("VCF [{}] was not in the set of recognized callers."
                          .format(vcf.name))

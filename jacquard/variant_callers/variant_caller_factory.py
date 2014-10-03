from jacquard.variant_callers.varscan import Varscan
from jacquard.variant_callers.strelka import Strelka
from jacquard.variant_callers.mutect import Mutect

from jacquard.jacquard_utils import JQException

_CALLERS = [Varscan(), Strelka(), Mutect()]

def get_caller(vcf):
    for caller in _CALLERS:
        #TODO: this should pass vcf instead of header
        if caller.validate_input_file(vcf.header):
            return caller
    raise JQException("VCF [{}] was not in the set of recognized callers."
                      .format(vcf.name))


#TODO cgates: suspect this does not need to be a class
class Factory(object):
    def __init__(self):
        self._callers = _CALLERS

    def get_vcf_caller(self, vcf):
        for caller in self._callers:
            #TODO: this should pass vcf instead of header
            if caller.validate_input_file(vcf.header):
                return caller
        raise JQException("VCF [{}] was not in the set of recognized callers."
                          .format(vcf.name))

    def get_caller(self, header):
        for caller in self._callers:
            if caller.validate_input_file(header):
                return caller
        raise JQCallerNotRecognizedException()

#TODO cgates: remove when get_caller is removed
class JQCallerNotRecognizedException(JQException):
    pass


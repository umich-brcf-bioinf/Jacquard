from jacquard.variant_callers.varscan import Varscan
from jacquard.variant_callers.strelka import Strelka
from jacquard.variant_callers.mutect import Mutect

import jacquard.utils as utils

_CALLERS = [Varscan(), Strelka(), Mutect()]

def get_caller(metaheaders, column_header, name):
    for caller in _CALLERS:
        #TODO: this should pass vcf instead of header
        if caller.validate_input_file(metaheaders, column_header):
            utils.log("DEBUG: VCF [{}] recognized by caller [{}]",
                 name, caller.name)
            return caller
    raise utils.JQException("VCF [{}] was not in the set of recognized callers."
                      .format(name))

@property
def callers():
    return _CALLERS

# def recognize_vcf(vcf_reader,):


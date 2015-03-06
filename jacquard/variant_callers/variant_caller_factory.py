"""Delegates to third-party variant callers.

To keep things simple, fair, and encapsulated route any requests to variant
callers through this module; likewise, don't directly reference individual
variant callers outside of this module.
"""
from __future__ import print_function, absolute_import
from jacquard.variant_callers.mutect import Mutect
from jacquard.variant_callers.strelka import Strelka
from jacquard.variant_callers.varscan import Varscan
import jacquard.logger as logger
import jacquard.utils as utils

#TODO: cgates: These should be defined by the caller modules themselves.
SUPPORTED_CALLER_VERSIONS = {"VarScan":"v2.3",
                             "MuTect": "v1.1.4",
                             "Strelka": "v2.0.15"}
_CALLERS = [Varscan(), Strelka(), Mutect()]

#TODO: (cgates): Filter uses this, but only for logging; adjust filter and drop
# method. Then consider renaming the module.
def get_caller(metaheaders, column_header, name):
    for caller in _CALLERS:
        if caller.validate_input_file(metaheaders, column_header):
            logger.debug("VCF [{}] recognized by caller [{}]",
                         name,
                         caller.name)
            return caller
    raise utils.JQException(("VCF [{}] was not in the set of "
                             "recognized callers.").format(name))

def claim(unclaimed_file_readers):
    """Allows each caller to claim incoming files as they are recognized.
    
    Args:
        unclaimed_file_readers: Usually, all files in the input dir.

    Returns:
        A tuple of unclaimed file readers and claimed VcfReaders. The presence
        of any unclaimed file readers could indicate stray files in the input
        dir.
    """
    claimed_vcf_readers = []
    for caller in _CALLERS:
        (unclaimed_file_readers,
         translated_vcf_readers) = caller.claim(unclaimed_file_readers)
        claimed_vcf_readers.extend(translated_vcf_readers)
    return unclaimed_file_readers, claimed_vcf_readers


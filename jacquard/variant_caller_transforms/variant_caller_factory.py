"""Delegates to third-party variant callers.

To keep things simple, fair, and encapsulated route any requests to variant
callers through this module; likewise, don't directly reference individual
variant callers outside of this module.
"""
from __future__ import print_function, absolute_import, division

import jacquard.utils.logger as logger
import jacquard.utils.utils as utils
import jacquard.variant_caller_transforms.mutect as mutect
import jacquard.variant_caller_transforms.strelka as strelka
import jacquard.variant_caller_transforms.varscan as varscan


SUPPORTED_CALLER_VERSIONS = {"VarScan": varscan.VERSION,
                             "MuTect": mutect.VERSION,
                             "Strelka": strelka.VERSION}

class VariantCallerFactory(object):
    def __init__(self, args=None):
        self._callers = [varscan.Varscan(args),
                         strelka.Strelka(),
                         mutect.Mutect()]

#TODO: (cgates): Filter uses this, but only for logging; adjust filter and drop
# method. Then consider renaming the module.
    def get_caller(self, metaheaders, column_header, name):
        for caller in self._callers:
            if caller.validate_input_file(metaheaders, column_header):
                logger.debug("VCF [{}] recognized by caller [{}]",
                             name,
                             caller.name)
                return caller
        raise utils.JQException(("VCF [{}] was not in the set of "
                                 "recognized callers.").format(name))

    def claim(self, unclaimed_file_readers):
        """Allows each caller to claim incoming files as they are recognized.

        Args:
            unclaimed_file_readers: Usually, all files in the input dir.

        Returns:
            A tuple of unclaimed file readers and claimed VcfReaders. The
            presence of any unclaimed file readers could indicate stray files
            in the input dir.
        """
        claimed_vcf_readers = []
        for caller in self._callers:
            (unclaimed_file_readers,
             translated_vcf_readers) = caller.claim(unclaimed_file_readers)
            claimed_vcf_readers.extend(translated_vcf_readers)

        return unclaimed_file_readers, claimed_vcf_readers

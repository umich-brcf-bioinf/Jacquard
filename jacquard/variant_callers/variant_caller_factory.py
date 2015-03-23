"""Delegates to third-party variant callers.

To keep things simple, fair, and encapsulated route any requests to variant
callers through this module; likewise, don't directly reference individual
variant callers outside of this module.
"""
from __future__ import print_function, absolute_import

import jacquard.logger as logger
import jacquard.utils as utils
import jacquard.variant_callers.mutect as mutect
import jacquard.variant_callers.strelka as strelka
import jacquard.variant_callers.varscan as varscan


SUPPORTED_CALLER_VERSIONS = {"VarScan": varscan.VERSION,
                             "MuTect": mutect.VERSION,
                             "Strelka": strelka.VERSION}

class VariantCallerFactory(object):
    def __init__(self, args=None):
        self.hc_filter_filename = None
        self.allow_inconsistent_sample_sets = None
        if args:
            self._handle_arguments(args)
        self._callers = [varscan.Varscan(args),
                         strelka.Strelka(),
                         mutect.Mutect()]

    def _handle_arguments(self, args):
        try:
            if args.varscan_hc_filter_filename:
                self.hc_filter_filename = args.varscan_hc_filter_filename
            if args.allow_inconsistent_sample_sets:
                #pylint: disable=line-too-long
                self.allow_inconsistent_sample_sets = args.allow_inconsistent_sample_sets
        except AttributeError:
            pass

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

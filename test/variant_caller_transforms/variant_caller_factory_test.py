#pylint:disable=line-too-long,invalid-name,global-statement, too-many-public-methods
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import re

import jacquard.variant_caller_transforms.variant_caller_factory as variant_caller_factory
import test.utils.test_case as test_case
import test.utils.vcf_test as vcf_test


class VariantCallerFactoryClassTestCase(test_case.JacquardBaseTestCase):
    def test_claim(self):
        args = Namespace(varscan_hc_filter_filename=None)
        factory = variant_caller_factory.VariantCallerFactory(args)
        fileA = vcf_test.MockFileReader("fileA.vcf")
        fileB = vcf_test.MockFileReader("fileB.vcf")
        fileC = vcf_test.MockFileReader("fileC.vcf")
        file_readers = [fileA,
                        fileB,
                        fileC]
        factory._callers = [vcf_test.MockCaller("foo", claimable=[fileC]),
                            vcf_test.MockCaller("bar", claimable=[fileB]),
                            vcf_test.MockCaller("baz", claimable=[fileA])]
        expected = [fileC, fileB, fileA]
        _, actual = factory.claim(file_readers)
        self.assertEquals(expected, actual)

    def test_claim_unclaimedFilesRemain(self):
        args = Namespace(varscan_hc_filter_filename=None)
        factory = variant_caller_factory.VariantCallerFactory(args)
        fileA = vcf_test.MockFileReader("fileA.vcf")
        fileB = vcf_test.MockFileReader("fileB.vcf")
        unclaimable_file = vcf_test.MockFileReader("fileC.txt")
        file_readers = [fileA,
                        fileB,
                        unclaimable_file]
        factory._callers = [vcf_test.MockCaller("foo", claimable=[fileB]),
                                           vcf_test.MockCaller("bar", claimable=[fileA])]
        expected_unclaimed = [unclaimable_file]
        actual_unclaimed, _ = factory.claim(file_readers)
        self.assertEquals(expected_unclaimed, actual_unclaimed)

    def test_claim_initializingVarscanStoresHCFile(self):
        args = Namespace(varscan_hc_filter_filename="bar$")
        factory = variant_caller_factory.VariantCallerFactory(args)
        self.assertEquals(re.compile("bar$"), factory._callers[0].hc_file_pattern)

#pylint:disable=line-too-long,invalid-name,global-statement, too-many-public-methods
import test.test_case as test_case
import jacquard.variant_callers.variant_caller_factory as variant_caller_factory
from jacquard.utils import JQException
import test.vcf_test as vcf_test

original_callers = None

def mock_callers():
    return [vcf_test.MockCaller("foo"),
            vcf_test.MockCaller("bar"),
            vcf_test.MockCaller("baz")]

class VariantCallerFactoryTestCase(test_case.JacquardBaseTestCase):
    def test_defined_caller(self):
        self.assertEquals(3, len(variant_caller_factory._CALLERS))
        caller_names = [caller.name for caller in variant_caller_factory._CALLERS]
        self.assertIn('VarScan', caller_names)
        self.assertIn('MuTect', caller_names)
        self.assertIn('Strelka', caller_names)

    def test_caller_notFoundRaisesException(self):
        self.assertRaises(JQException,
                          variant_caller_factory.get_caller,
                          ["##metaheaders"],
                          "#header",
                          "vcfName")

class VariantCallerFactoryClaimTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(VariantCallerFactoryClaimTestCase, self).setUp()
        global original_callers
        original_callers = list(variant_caller_factory._CALLERS)
        variant_caller_factory._CALLERS = mock_callers()

    def tearDown(self):
        variant_caller_factory._CALLERS = original_callers
        super(VariantCallerFactoryClaimTestCase, self).tearDown()

    def test_claim(self):
        file_readers = [vcf_test.MockFileReader("fileA.vcf"),
                        vcf_test.MockFileReader("fileB.vcf"),
                        vcf_test.MockFileReader("fileC.vcf")]
        expected = ["foo", "bar", "baz"]
        actual = variant_caller_factory.claim(file_readers)
        self.assertEquals(expected, actual)

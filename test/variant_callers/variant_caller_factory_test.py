#pylint: disable=missing-docstring,line-too-long,too-many-public-methods
#pylint: disable=too-few-public-methods,too-many-instance-attributes
#pylint: disable=too-many-arguments,invalid-name,protected-access,global-statement
import unittest
import jacquard.variant_callers.variant_caller_factory as variant_caller_factory
from jacquard.utils import JQException


class VariantCallerFactoryTestCase(unittest.TestCase):
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

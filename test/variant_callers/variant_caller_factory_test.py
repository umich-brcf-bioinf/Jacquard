import unittest
import jacquard.variant_callers.variant_caller_factory as variant_caller_factory
from jacquard.jacquard_utils import JQException

class VariantCallerFactoryTestCase(unittest.TestCase):
    def test_defined_caller(self):
        self.assertEquals(3, len(variant_caller_factory._CALLERS))
        caller_names = [caller.name for caller in variant_caller_factory._CALLERS]
        self.assertIn('VarScan', caller_names)
        self.assertIn('MuTect', caller_names)
        self.assertIn('Strelka', caller_names)
    
    def test_caller_notFoundRaisesException(self):
        self.assertRaises(JQException, variant_caller_factory.get_caller, ["##unknown metaheaders"], "#column_header", "foo")
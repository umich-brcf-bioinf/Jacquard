import unittest

from jacquard.variant_callers import unknown

class UnknownTestCase(unittest.TestCase):
    def test_validateInputFile_valid(self):
        caller = unknown.Unknown() 
        input_file = ["##bar", "foo"]
        name, valid = caller.validate_input_file(input_file)
        
        self.assertEquals("Unknown", name)
        self.assertEquals(1, valid)
        

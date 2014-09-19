from collections import OrderedDict,defaultdict
import glob
import os
from os import listdir
from StringIO import StringIO
import sys
from stat import *
import testfixtures
from testfixtures import TempDirectory
import unittest

import unknown

class UnknownTestCase(unittest.TestCase):
    def test_validateInputFile_valid(self):
        caller = unknown.Unknown() 
        input_file = ["##bar", "foo"]
        name, valid = caller.validate_input_file(input_file)
        
        self.assertEquals("Unknown", name)
        self.assertEquals(1, valid)
        

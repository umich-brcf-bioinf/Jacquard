'''
Created on Oct 22, 2014

@author: kmeng
'''
import unittest
import jacquard.logger as logger
from testfixtures import TempDirectory
import sys
from StringIO import StringIO
from datetime import datetime

class LoggerTestCase(unittest.TestCase):


    def setUp(self):
        self.output_dir = TempDirectory()
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output
        
    def tearDown(self):
        self.output_dir.cleanup()
        self.output.close()
        sys.stderr = self.saved_stderr

    def test_initialize_logger(self):
        tool = "foo"
        logger.initialize_logger(self.output_dir.path,tool)
        
        self.assertEquals(['host','tool','user','time'], logger.logging_dict.keys())
    
    def test_error(self):
        tool = "foo"
        logger.initialize_logger(self.output_dir.path,tool)
        logger.error("bar")
        current_time = datetime.now().strftime('%Y/%m/%d')
        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|ERROR\|foo\|bar")
         
    def test_warning(self):
        tool = "foo"
        logger.initialize_logger(self.output_dir.path,tool)
        logger.warning("bar")
        current_time = datetime.now().strftime('%Y/%m/%d')
        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|WARNING\|foo\|bar")
         
    def test_info(self):
        tool = "foo"
        logger.initialize_logger(self.output_dir.path,tool)
        logger.info("bar")
        current_time = datetime.now().strftime('%Y/%m/%d')
        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|INFO\|foo\|bar")
 
    def test_debug(self):
        tool = "foo"
        logger.initialize_logger(self.output_dir.path,tool)
        logger.debug("bar")
        current_time = datetime.now().strftime('%Y/%m/%d')
        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|DEBUG\|foo\|bar")
    
#     def test_validation_messages(self):
#         tool = "foo"
#         logger.initialize_logger(self.output_dir.path,tool)
#         logger.validation_messages("bar","baz")
#         current_time = datetime.now().strftime('%Y/%m/%d')
#         output_lines = self.output.getvalue().rstrip().split("\n")
#         self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|WARNING\|foo\|baz")
#         self.assertRegexpMatches(output_lines[1], ""+current_time+r".*\|ERROR\|foo\|bar")
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
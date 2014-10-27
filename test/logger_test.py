'''
Created on Oct 22, 2014

@author: kmeng
'''
from datetime import datetime
import os
import sys
from StringIO import StringIO
import unittest
import jacquard.logger as logger 

class LoggerTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
#         self.saved_stderr = sys.stderr
        sys.stderr = self.output
        self.log_file = os.path.join(os.path.dirname(os.getcwd()), "logs", "jacquard.log")
        try:
            os.remove(self.log_file)
        except:
            pass
        
    def tearDown(self):
        self.output.close()
#         sys.stderr = self.saved_stderr
#         try:
#             os.remove(self.log_file)
#         except:
#             pass

    def test_initialize_logger(self):
        tool = "foo"
        logger.initialize_logger(tool)
        
        self.assertEquals(['host','tool','user','time'], logger.logging_dict.keys())
    
    def test_error(self): 
        tool = "foo"
        logger.initialize_logger(tool)
        logger.error("bar")
        root_logger = logger.logging.getLogger()
        
        current_time = datetime.now().strftime('%Y/%m/%d')
        output_lines = self.output.getvalue().rstrip().split("\n")

        self.assertEquals(["root: ERROR: bar"], root_logger.handlers[0].buffer)
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|ERROR\|foo\|bar")
         
    def test_warning(self):
        tool = "foo"
        logger.initialize_logger(tool)
        logger.warning("bar")
        root_logger = logger.logging.getLogger()
        
        current_time = datetime.now().strftime('%Y/%m/%d')
        output_lines = self.output.getvalue().rstrip().split("\n")
        
        self.assertEquals(["root: WARNING: bar"], root_logger.handlers[0].buffer)
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|WARNING\|foo\|bar")
         
    def test_info(self):
        tool = "foo"
        logger.initialize_logger(tool)
        logger.info("bar")
        root_logger = logger.logging.getLogger()
        
        current_time = datetime.now().strftime('%Y/%m/%d')
        output_lines = self.output.getvalue().rstrip().split("\n")
        
        self.assertEquals(["root: INFO: bar"], root_logger.handlers[0].buffer)
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|INFO\|foo\|bar")
 
    def test_debug(self):
        tool = "foo"
        logger.initialize_logger(tool)
        logger.debug("bar")
        root_logger = logger.logging.getLogger()
        
        current_time = datetime.now().strftime('%Y/%m/%d')
        output_lines = self.output.getvalue().rstrip().split("\n")
#         self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|DEBUG\|foo\|bar")
    
        self.assertEquals(["root: DEBUG: bar"], root_logger.handlers[0].buffer)
        self.assertEquals(output_lines, [""])
    
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
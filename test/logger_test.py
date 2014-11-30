# pylint: disable=C0103,C0301,R0903,R0904
from datetime import datetime
import os
import sys
from StringIO import StringIO
import unittest
import jacquard.logger as logger


class LoggerTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        sys.stderr = self.output

    def tearDown(self):
        self.output.close()
        log_file = os.path.join(os.path.dirname(os.getcwd()), "jacquard.log")
        try:
            os.remove(log_file)
        except OSError:
            pass

    def test_initialize_logger(self):
        tool = "foo"
        logger.initialize_logger(tool)
        self.assertEquals(['host', 'tool', 'start_time', 'user'],
                          logger.logging_dict.keys())

    def test_error(self):
        tool = "foo"
        logger.initialize_logger(tool)
        logger.error("bar")

        root_logger = logger.logging.getLogger()

        current_time = datetime.now().strftime('%Y-%m-%d')
        output_lines = self.output.getvalue().rstrip().split("\n")

        ##nosetests overwrites logger.FileHandler
        self.assertEquals(["root: ERROR: bar"], root_logger.handlers[0].buffer)

        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|ERROR\|foo\|bar")

    def test_warning(self):
        tool = "foo"
        logger.initialize_logger(tool)
        self.assertFalse(logger.SHOW_WARNING)
        logger.warning("bar")
        self.assertTrue(logger.SHOW_WARNING)
        root_logger = logger.logging.getLogger()

        current_time = datetime.now().strftime('%Y-%m-%d')
        output_lines = self.output.getvalue().rstrip().split("\n")

        self.assertEquals(["root: WARNING: bar"], root_logger.handlers[0].buffer)
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|WARNING\|foo\|bar")

    def test_info(self):
        tool = "foo"
        logger.initialize_logger(tool)
        current_time = datetime.now()

        logger.info("bar")

        root_logger = logger.logging.getLogger()
        self.assertEquals(["root: INFO: bar"], root_logger.handlers[0].buffer)

        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertEquals(len(output_lines), 1)

        log_fields = output_lines[0].split("|")
        self.assertEquals(4, len(log_fields))
        log_time_string = log_fields[0]
        log_datetime = datetime.strptime(log_time_string, "%Y-%m-%d %H:%M:%S")
        seconds_delta = (log_datetime - current_time).total_seconds()
        self.assertLess(seconds_delta, 1)
        self.assertEquals(["INFO", "foo", "bar"], log_fields[1:])

    def test_info_message_args(self):
        tool = "foo"
        logger.initialize_logger(tool)
        current_time = datetime.now()

        logger.info("bar {}/{}", "1", "2")

        root_logger = logger.logging.getLogger()
        self.assertEquals(["root: INFO: bar 1/2"], root_logger.handlers[0].buffer)

        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertEquals(len(output_lines), 1)

        log_fields = output_lines[0].split("|")
        self.assertEquals(4, len(log_fields))
        log_time_string = log_fields[0]
        log_datetime = datetime.strptime(log_time_string, "%Y-%m-%d %H:%M:%S")
        seconds_delta = (log_datetime - current_time).total_seconds()
        self.assertLess(seconds_delta, 1)
        self.assertEquals(["INFO", "foo", "bar 1/2"], log_fields[1:])

    def test_debug(self):
        tool = "foo"
        logger.initialize_logger(tool)

        logger.debug("bar")

        root_logger = logger.logging.getLogger()
        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertEquals(["root: DEBUG: bar"], root_logger.handlers[0].buffer)
        self.assertEquals(output_lines, [""])

    def test_debug_consoleWhenVerbose(self):
        tool = "foo"
        logger.initialize_logger(tool, True)
        logger.debug("bar")
        root_logger = logger.logging.getLogger()

        current_time = datetime.now().strftime('%Y-%m-%d')
        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|DEBUG\|foo\|bar")

        self.assertEquals(["root: DEBUG: bar"], root_logger.handlers[0].buffer)
        self.assertRegexpMatches(output_lines[0], r""+current_time+"|DEBUG|foo|bar'")

    def test_debug_noConsoleWhenNotVerbose(self):
        tool = "foo"
        logger.initialize_logger(tool)

        logger.debug("bar")

        output_lines = self.output.getvalue().rstrip().split("\n")
        root_logger = logger.logging.getLogger()
        self.assertEquals(["root: DEBUG: bar"], root_logger.handlers[0].buffer)
        self.assertEquals(output_lines[0], "")

    def test_noExceptionOnMalformedMessage(self):
        tool = "foo"
        logger.initialize_logger(tool)

        logger.info("bar {}/{}/{}", "1", "2")

        root_logger = logger.logging.getLogger()
        self.assertEquals(["root: INFO: Malformed log message (IndexError: tuple index out of range)|bar {}/{}/{}|['1', '2']"], root_logger.handlers[0].buffer)

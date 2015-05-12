#pylint: disable=line-too-long, invalid-name, too-many-public-methods
from __future__ import print_function, absolute_import, division

from argparse import Namespace
from datetime import datetime
import os
import shutil
import sys
import unittest

import jacquard.utils.logger as logger

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


class LoggerTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        sys.stderr = self.output
        logger._verbose = False

    def tearDown(self):
        self.output.close()
        log_file = os.path.join(os.path.dirname(os.getcwd()), "jacquard.log")
        try:
            os.remove(log_file)
        except OSError:
            pass
        logger._verbose = False

    def test_initialize_logger_defaultFilename(self):
        tool = "foo"
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=None)
        logger.initialize_logger(args)
        self.assertEquals(['host', 'start_time', 'tool', 'user'],
                          sorted(logger._logging_dict.keys()))
        self.assertEquals("jacquard.log", os.path.basename(logger.log_filename))

    def test_initialize_logger_suppliedFilename(self):
        try:
            tool = "foo"
            log_filename = "tmp/log.foo"
            args = Namespace(subparser_name=tool,
                             log_file=log_filename,
                             verbose=None)
            logger.initialize_logger(args)
            self.assertEquals(['host', 'start_time', 'tool', 'user'],
                              sorted(logger._logging_dict.keys()))
            self.assertEquals(log_filename, logger.log_filename)

        finally:
            shutil.rmtree(os.path.dirname(log_filename))

    def test_initialize_logger_createsParentDirs(self):
        try:
            tool = "foo"
            log_filename = "tmp/log.foo"
            args = Namespace(subparser_name=tool,
                             log_file=log_filename,
                         verbose=None)
            logger.initialize_logger(args)
            self.assertTrue(os.path.isdir(os.path.dirname(log_filename)))
        finally:
            shutil.rmtree(os.path.dirname(log_filename))

    def test_initialize_logger_verbose(self):
        tool = "foo"
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=True)
        logger.initialize_logger(args)
        logger.debug("bar")
        root_logger = logger.logging.getLogger()

        current_time = datetime.now().strftime('%Y-%m-%d')
        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|DEBUG\|foo\|bar")

        self.assertEquals(["root: DEBUG: bar"], root_logger.handlers[0].buffer)
        self.assertRegexpMatches(output_lines[0], r""+current_time+"|DEBUG|foo|bar'")

    def test_initialize_logger_notVerbose(self):
        tool = "foo"
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=None)
        logger.initialize_logger(args)
        logger.debug("bar")

        output_lines = self.output.getvalue().rstrip().split("\n")
        root_logger = logger.logging.getLogger()
        self.assertEquals(["root: DEBUG: bar"], root_logger.handlers[0].buffer)
        self.assertEquals(output_lines[0], "")

    def test_error(self):
        tool = "foo"
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=None)
        logger.initialize_logger(args)
        logger.error("bar")

        root_logger = logger.logging.getLogger()

        current_time = datetime.now().strftime('%Y-%m-%d')
        output_lines = self.output.getvalue().rstrip().split("\n")

        ##nosetests overwrites logger.FileHandler
        self.assertEquals(["root: ERROR: bar"], root_logger.handlers[0].buffer)

        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|ERROR\|foo\|bar")

    def test_warning(self):
        tool = "foo"
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=None)
        logger.initialize_logger(args)
        self.assertFalse(logger.WARNING_OCCURRED)
        logger.warning("bar")
        self.assertTrue(logger.WARNING_OCCURRED)
        root_logger = logger.logging.getLogger()

        current_time = datetime.now().strftime('%Y-%m-%d')
        output_lines = self.output.getvalue().rstrip().split("\n")

        self.assertEquals(["root: WARNING: bar"], root_logger.handlers[0].buffer)
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|WARNING\|foo\|bar")

    def test_info(self):
        tool = "foo"
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=None)
        logger.initialize_logger(args)
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
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=None)
        logger.initialize_logger(args)
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

    def test_noExceptionOnMalformedMessage(self):
        tool = "foo"
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=None)
        logger.initialize_logger(args)

        logger.info("bar {}/{}/{}", "1", "2")

        root_logger = logger.logging.getLogger()
        self.assertEquals(["root: INFO: Malformed log message (IndexError: tuple index out of range)|bar {}/{}/{}|['1', '2']"], root_logger.handlers[0].buffer)

class LoggerTestCaseDebugOnly(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        sys.stderr = self.output
        logger._verbose = False

    def tearDown(self):
        self.output.close()
        log_file = os.path.join(os.path.dirname(os.getcwd()), "jacquard.log")
        try:
            os.remove(log_file)
        except OSError:
            pass
        logger._verbose = False

    def test_debug(self):
        tool = "foo"
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=None)
        logger.initialize_logger(args)

        logger.debug("bar")

        root_logger = logger.logging.getLogger()
        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertEquals(["root: DEBUG: bar"], root_logger.handlers[0].buffer)
        self.assertEquals(output_lines, [""])

    def test_debug_noConsoleWhenNotVerbose(self):
        tool = "foo"
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=None)
        logger.initialize_logger(args)

        logger.debug("bar")

        output_lines = self.output.getvalue().rstrip().split("\n")
        root_logger = logger.logging.getLogger()
        self.assertEquals(["root: DEBUG: bar"], root_logger.handlers[0].buffer)
        self.assertEquals(output_lines[0], "")

    def test_debug_consoleWhenVerbose(self):
        tool = "foo"
        args = Namespace(subparser_name=tool,
                         log_file=None,
                         verbose=True)
        logger.initialize_logger(args)
        logger.debug("bar")
        root_logger = logger.logging.getLogger()

        current_time = datetime.now().strftime('%Y-%m-%d')
        output_lines = self.output.getvalue().rstrip().split("\n")
        self.assertRegexpMatches(output_lines[0], ""+current_time+r".*\|DEBUG\|foo\|bar")

        self.assertEquals(["root: DEBUG: bar"], root_logger.handlers[0].buffer)
        self.assertRegexpMatches(output_lines[0], r""+current_time+"|DEBUG|foo|bar'")

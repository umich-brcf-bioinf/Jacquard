#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import os
import shutil
import sys
import unittest

import jacquard.jacquard as jacquard
import jacquard.utils.logger as logger
import jacquard.utils.vcf as vcf


try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO



class JacquardBaseTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output

    def tearDown(self):
        logger.WARNING_OCCURRED = False
        self.output.close()
        sys.stderr = self.saved_stderr
        unittest.TestCase.tearDown(self)

    def assertStartsWith(self, full_text, search_text):
        self.assertTrue(full_text.startswith(search_text))

    def ok(self):
        self.assertTrue(True)

    def move_files(self, source_dirs, dest_dir):
        try:
            os.mkdir(dest_dir)
        except IOError:
            pass
        for source_dir in source_dirs:
            for output_file in os.listdir(source_dir):
                shutil.move(os.path.join(source_dir, output_file),
                            dest_dir)

    def _compare_files(self, actual_output_path, expected_output_path):
        actual_file = vcf.FileReader(actual_output_path)
        actual_file.open()
        actual = [line for line in actual_file.read_lines()]
        actual_file.close()

        expected_file = vcf.FileReader(expected_output_path)
        expected_file.open()
        expected = [line for line in expected_file.read_lines()]
        expected_file.close()

        self.assertEquals(len(expected),
                          len(actual))

        for i in range(len(expected)):
            if expected[i].startswith("##jacquard=<timestamp="):
                self.assertStartsWith(actual[i], "##jacquard=<timestamp=")
            else:
                self.assertEquals(expected[i].rstrip(),
                                  actual[i].rstrip())

    def assertCommand(self, command, expected_output_path):
        jacquard._dispatch(jacquard._SUBCOMMANDS, command)

        actual_output_path = command[2]

        try:
            if os.path.isfile(actual_output_path):
                self._compare_files(actual_output_path, expected_output_path)
            elif os.path.isdir(actual_output_path):
                for output_filename in os.listdir(actual_output_path):
                    actual_output_file = os.path.join(actual_output_path,
                                                      output_filename)
                    expected_output_file = os.path.join(expected_output_path,
                                                        output_filename)
                    self._compare_files(actual_output_file,
                                        expected_output_file)
        except self.failureException as e:
            msg = "discrepancy in command [{}]: {}".format(" ".join(command), e)
            raise self.failureException(msg)

    def entab(self, string, old="|"):
        return string.replace(old, "\t")

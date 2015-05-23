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

    def _compare_files(self, output, output_file, expected_dir):
        actual_file = vcf.FileReader(output)
        actual_file.open()
        actual = [line for line in actual_file.read_lines()]
        actual_file.close()

        expected_file = vcf.FileReader(os.path.join(expected_dir,
                                                    output_file))
        expected_file.open()
        expected = [line for line in expected_file.read_lines()]
        expected_file.close()

        self.assertEquals(len(expected),
                          len(actual))

        for i in range(len(expected)):
            if expected[i].startswith("##jacquard=<Timestamp="):
                self.assertStartsWith(actual[i], "##jacquard=<Timestamp=")
            else:
                self.assertEquals(expected[i].rstrip(),
                                  actual[i].rstrip())

    def assertCommand(self, command, expected_dir):
        jacquard._dispatch(jacquard._SUBCOMMANDS, command)

        output = command[2]

        try:
            if os.path.isfile(output):
                output_file = os.path.basename(output)
                self._compare_files(output, output_file, expected_dir)
            elif os.path.isdir(output):
                for output_file in os.listdir(output):
                    new_output = os.path.join(output, output_file)
                    self._compare_files(new_output, output_file, expected_dir)
        except self.failureException as e:
            msg = "discrepancy in command [{}]: {}".format(" ".join(command), e)
            raise self.failureException(msg)

    def entab(self, string, old="|"):
        return string.replace(old, "\t")

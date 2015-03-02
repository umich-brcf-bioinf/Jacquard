#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import
import os
import shutil
from StringIO import StringIO
import sys
import unittest

import jacquard.jacquard as jacquard
import jacquard.logger as logger
import jacquard.vcf as vcf

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

        print(actual)

        self.assertEquals(len(expected),
                          len(actual))

        for i in xrange(len(expected)):
            if expected[i].startswith("##jacquard.cwd="):
                self.assertStartsWith(actual[i], "##jacquard.cwd=")
            elif expected[i].startswith("##jacquard.command="):
                self.assertStartsWith(actual[i], "##jacquard.command=")
            else:
                self.assertEquals(expected[i].rstrip(),
                                  actual[i].rstrip())

    def assertCommand(self, command, expected_dir):
        jacquard.dispatch(jacquard._SUBCOMMANDS, command)

        output = command[2]

        try:
            if os.path.isfile(output):
                output_file = os.path.basename(output)
                self._compare_files(output, output_file, expected_dir)
            elif os.path.isdir(output):
                for output_file in os.listdir(output):
                    new_output = os.path.join(output, output_file)
                    self._compare_files(new_output, output_file, expected_dir)
        except AssertionError as e:
            msg = "discrepancy in command [{}]: {}".format(" ".join(command), e)
            raise AssertionError(msg)

    def entab(self, string, old="|"):
        return string.replace(old, "\t")

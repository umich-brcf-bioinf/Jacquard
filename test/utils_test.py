#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=missing-docstring, protected-access, global-statement, too-few-public-methods
from __future__ import absolute_import, print_function

from StringIO import StringIO
import os
import subprocess
import sys

import natsort
from testfixtures import TempDirectory

import jacquard.logger as logger
import jacquard.utils as utils
import test.test_case as test_case

mock_log_called = False
mock_message = ""

#pylint: disable=unused-argument
def mock_log(msg, *args):
    global mock_log_called
    mock_log_called = True

class JQExceptionTestCase(test_case.JacquardBaseTestCase):
    def test_init(self):
        actual = utils.JQException("msg:{}, {}", "bar", [1, 2, 3])
        self.assertIsInstance(actual, Exception)
        self.assertEquals(actual.message, "msg:bar, [1, 2, 3]")

class OrderedSetTestCase(test_case.JacquardBaseTestCase):
    def test_isOrderedSet(self):
        actual = utils.OrderedSet(["B", "A", "B"])
        self.assertEquals(2, len(actual))
        it = iter(actual)
        self.assertEquals("B", it.next())
        self.assertEquals("A", it.next())

    def test_contains(self):
        actual = utils.OrderedSet(["A", "B"])
        self.assertIn("A", actual)
        self.assertIn("B", actual)
        self.assertNotIn("C", actual)

    def test_add(self):
        actual = utils.OrderedSet(["B", "A"])
        actual.add("A")
        self.assertEquals(2, len(actual))
        it = iter(actual)
        self.assertEquals("B", it.next())
        self.assertEquals("A", it.next())

        actual.add("C")
        self.assertEquals(3, len(actual))
        it = iter(actual)
        self.assertEquals("B", it.next())
        self.assertEquals("A", it.next())
        self.assertEquals("C", it.next())

    def test_discard(self):
        actual = utils.OrderedSet(["B", "A"])
        actual.discard("B")
        self.assertEquals(1, len(actual))
        it = iter(actual)
        self.assertEquals("A", it.next())

        actual.add("B")
        self.assertEquals(2, len(actual))
        it = iter(actual)
        self.assertEquals("A", it.next())
        self.assertEquals("B", it.next())

        actual.discard("C")
        self.assertEquals(2, len(actual))

    def test_reversed(self):
        actual = utils.OrderedSet(["B", "A"])

        self.assertEquals(2, len(actual))
        it = reversed(actual)
        self.assertEquals("A", it.next())
        self.assertEquals("B", it.next())

    def test_pop_right(self):
        actual = utils.OrderedSet(["B", "A"])

        self.assertEquals(2, len(actual))
        self.assertEquals("A", actual.pop())
        self.assertEquals(1, len(actual))
        self.assertEquals("B", actual.pop())
        self.assertEquals(0, len(actual))
        self.assertRaises(KeyError, actual.pop)

    def test_pop_left(self):
        actual = utils.OrderedSet(["B", "A"])

        self.assertEquals(2, len(actual))
        self.assertEquals("B", actual.pop(last=False))
        self.assertEquals(1, len(actual))
        self.assertEquals("A", actual.pop(last=False))
        self.assertEquals(0, len(actual))
        self.assertRaises(KeyError, actual.pop)

    def test_eq(self):
        base = utils.OrderedSet(["A", "B"])
        same = utils.OrderedSet(["A", "B"])
        equivalentSet = set(["A", "B"])
        equivalentList = ["A", "B"]
        differentClass = "foo"
        differentOrder = utils.OrderedSet(["B", "A"])
        differentMembers = utils.OrderedSet(["A", "B", "C"])

        self.assertEquals(base, same)
        self.assertEquals(base, equivalentSet)
        self.assertEquals(base, equivalentList)
        self.assertNotEquals(base, differentClass)
        self.assertNotEquals(base, differentOrder)
        self.assertNotEquals(base, differentMembers)

    def test_repr(self):
        actual = utils.OrderedSet([])
        self.assertRegexpMatches(actual.__repr__(), r"OrderedSet")
        actual.add("B")
        actual.add("A")
        self.assertRegexpMatches(actual.__repr__(), r"['B','A']")

class NaturalSortTestCase(test_case.JacquardBaseTestCase):
    def test_natsort(self):
        unsorted = ["123a", "1abc", "13d"]
        expected = ["1abc", "13d", "123a"]
        actual = natsort.natsorted(unsorted)
        self.assertEquals(expected, actual)

    def test_natsort_lowerAndUpperCase(self):
        unsorted = ["123ABC", "123abc", "1abc", "13d"]
        expected = ["1abc", "13d", "123ABC", "123abc"]
        actual = natsort.natsorted(unsorted)
        self.assertEquals(expected, actual)

    def test_natsort_baseAlphaSort(self):
        unsorted = ["A100", "B1", "C10", "D"]
        expected = ["A100", "B1", "C10", "D"]
        actual = natsort.natsorted(unsorted)
        self.assertEquals(expected, actual)

    def test_natsort_numericOrder(self):
        unsorted = ["B100", "B1", "B10", "A101"]
        expected = ["A101", "B1", "B10", "B100"]
        actual = natsort.natsorted(unsorted)
        self.assertEquals(expected, actual)

    def test_natsort_breaksTiesByAlpha(self):
        unsorted = ["X100B", "X100C", "X100A", "X10"]
        expected = ["X10", "X100A", "X100B", "X100C"]
        actual = natsort.natsorted(unsorted)
        self.assertEquals(expected, actual)


class ValidateDirectoriesTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        self.original_info = logger.info
        self.original_error = logger.error
        self.original_warning = logger.warning
        self.original_debug = logger.debug
        self._change_mock_logger()
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output

    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr
        self._reset_mock_logger()

    @staticmethod
    def _change_mock_logger():
        global mock_log_called
        mock_log_called = False
        logger.info = mock_log
        logger.error = mock_log
        logger.warning = mock_log
        logger.debug = mock_log

    def _reset_mock_logger(self):
        logger.info = self.original_info
        logger.error = self.original_error
        logger.warning = self.original_warning
        logger.debug = self.original_debug

    def test_validateDirectories_inputDirectoryDoesntExist(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_file = script_dir + "/functional_tests/utils_test/tag_varscan_test/foo"
        output_file = script_dir + "/functional_tests/utils_test/tag_varscan_test/output"

        with self.assertRaises(SystemExit) as cm:
            utils.validate_directories(input_file, output_file)
        self.assertEqual(cm.exception.code, 1)

        self.assertTrue(mock_log_called)
#         self.assertRegexpMatches(self.output.getvalue(),
#                                  r"Specified input directory \[.*\] does not exist.")

    def test_validateDirectories_outputDirectoryNotCreated(self):
        with TempDirectory() as input_file, TempDirectory() as output_file:
            input_file.write("A.txt",
                            "##source=VarScan2\n#CHROM\tNORMAL\tTUMOR\n")
            unwriteable_dir = os.path.join(output_file.path, "unwriteable")
            desired_dir = os.path.join(unwriteable_dir, "bar")

            try:
                make_unwritable_dir(unwriteable_dir)
                with self.assertRaises(SystemExit) as cm:
                    utils.validate_directories(input_file.path, desired_dir)

            finally:
                cleanup_unwriteable_dir(unwriteable_dir)

            self.assertEqual(cm.exception.code, 1)
            self.assertTrue(mock_log_called)

#     def test_validate_arguments_valid(self):
#         output_files = ["foo.vcf", "bar.vcf"]
#         writers = [MockWriter(output_filepath="baz.vcf")]
#         self.assertTrue(utils.validate_arguments(output_files, writers))

#     def test_validate_arguments_invalid(self):
#         output_files = ["foo.vcf", "bar.vcf"]
#         writers = [MockWriter(output_filepath="bar.vcf"), MockWriter(output_filepath="blah.vcf")]
#         self.assertRaisesRegexp(utils.JQException, "This command would overwrite existing files", utils.validate_arguments, output_files, writers)

def is_windows_os():
    return sys.platform.lower().startswith("win")

def make_unwritable_dir(unwriteable_dir):
    os.mkdir(unwriteable_dir, 0555)

    if is_windows_os():
        FNULL = open(os.devnull, 'w')
        subprocess.call("icacls {0} /deny Everyone:W".format(unwriteable_dir), stdout=FNULL, stderr=subprocess.STDOUT)

def cleanup_unwriteable_dir(unwriteable_dir):
    if is_windows_os():
        FNULL = open(os.devnull, 'w')
        subprocess.call("icacls {0} /reset /t /c".format(unwriteable_dir), stdout=FNULL, stderr=subprocess.STDOUT)
    os.rmdir(unwriteable_dir)

class MockWriter(object):
    def __init__(self, output_filepath=None):
        self._content = []
        self.output_filepath = output_filepath
        self.wasClosed = False

    def write(self, content):
        self._content.extend(content.splitlines())

    def lines(self):
        return self._content

    def close(self):
        self.wasClosed = True


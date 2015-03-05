#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=missing-docstring, protected-access, global-statement, too-few-public-methods
from __future__ import absolute_import, print_function
import jacquard.utils as utils
import natsort
import test.test_case as test_case


class RoundTwoDigitsTestCase(test_case.JacquardBaseTestCase):
    def test_round_two_digits_noRounding(self):
        val = "1.01"
        expected = "1.01"
        actual = utils.round_two_digits(val)
        self.assertEquals(expected, actual)

        val = "1.1"
        expected = "1.1"
        actual = utils.round_two_digits(val)
        self.assertEquals(expected, actual)

    def test_round_two_digits_rounding(self):
        val = "1.011"
        expected = "1.01"
        actual = utils.round_two_digits(val)
        self.assertEquals(expected, actual)

        val = "1.016"
        expected = "1.02"
        actual = utils.round_two_digits(val)
        self.assertEquals(expected, actual)


class JQExceptionTestCase(test_case.JacquardBaseTestCase):
    def test_init(self):
        actual = utils.JQException("msg:{}, {}", "bar", [1, 2, 3])
        self.assertIsInstance(actual, Exception)
        self.assertEquals(actual.message, "msg:bar, [1, 2, 3]")


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


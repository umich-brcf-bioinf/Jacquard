#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=missing-docstring, protected-access, global-statement, too-few-public-methods
from __future__ import print_function, absolute_import, division

import argparse
import natsort

import jacquard.utils.utils as utils
import test.utils.test_case as test_case


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

class JaquardHelpFormatterTestCase(test_case.JacquardBaseTestCase):
    def test_extends_raw_text_help_formatter(self):
        jq_formatter = utils._JacquardHelpFormatter("prog")
        actual_lines = jq_formatter._split_lines("foo\nbar\nbaz", 100)
        expected_lines = ["foo", "bar", "baz"]

        self.assertEquals(expected_lines, actual_lines)

    def test_format_usage(self):
        jq_formatter = utils._JacquardHelpFormatter("prog")
        default_values = "[--include_rows=valid] [--foo=bar]"
        actual_usage = jq_formatter._format_usage(default_values)
        expected_usage = "usage: prog <input> <output> [--include_rows=valid] [--foo=bar]"

        self.assertEquals(expected_usage, actual_usage)

    def test_add_usage(self):
        jq_formatter = utils._JacquardHelpFormatter("prog")
        default_values = ["[--include_rows=valid]", "[--foo=bar]"]
        jq_formatter.add_usage(default_values)
        actual_tuple = jq_formatter._current_section.items

        self.assertEquals((jq_formatter._format_usage, default_values), actual_tuple[0])

    def test_jacquard_formatter(self):
        default_values = ["[--include_rows=valid]", "[--foo=bar]"]
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(title="subcommands",
                                           dest="subparser_name")
        parser = subparsers.add_parser("foo",
                                       formatter_class=utils._JacquardHelpFormatter,
                                       usage=default_values)
        self.assertEquals(default_values, parser.usage)

class JQExceptionTestCase(test_case.JacquardBaseTestCase):
    def test_init(self):
        actual = utils.JQException("msg:{}, {}", "bar", [1, 2, 3])
        self.assertIsInstance(actual, Exception)
        self.assertEquals(actual.args[0], "msg:bar, [1, 2, 3]")

class SortMetaheadersTestCase(test_case.JacquardBaseTestCase):
    def test_sort_metaheaders(self):
        unsorted = ["##FORMAT", "##fileformat", "##INFO", "##jacquard.foo="]
        actual = utils.sort_metaheaders(unsorted)
        expected = ["##fileformat", "##jacquard.foo=", "##INFO", "##FORMAT"]
        self.assertEquals(expected, actual)

    def test_sort_metaheaders_completeList(self):
        unsorted = ["##FORMAT", "##fileformat", "##INFO", "##jacquard.foo=", "##ALT", "#CHROM", "##contig", "##FILTER"]
        actual = utils.sort_metaheaders(unsorted)
        expected = ["##fileformat", "##jacquard.foo=", "##contig", "##ALT", "##FILTER", "##INFO", "##FORMAT", "#CHROM"]
        self.assertEquals(expected, actual)

    def test_sort_metaheaders_sortsWithinCategory(self):
        unsorted = ["##FORMAT=DP", "##FORMAT=AF", "##jacquard.bar="]
        actual = utils.sort_metaheaders(unsorted)
        expected = ["##jacquard.bar=", "##FORMAT=AF", "##FORMAT=DP"]
        self.assertEquals(expected, actual)

    def test_sort_metaheaders_sortsWithinCategoryAlphaNum(self):
        unsorted = ["##contig=<ID=chr13", "##contig=<ID=chr1", "##contig=<ID=chr2"]
        actual = utils.sort_metaheaders(unsorted)
        expected = ["##contig=<ID=chr1", "##contig=<ID=chr2", "##contig=<ID=chr13"]
        self.assertEquals(expected, actual)

    def test_sort_metaheaders_unexpectedMetaheadesr(self):
        unsorted = ["##FORMAT", "##INFO", "##Strelka", "##MuTect", "##VarScan"]
        actual = utils.sort_metaheaders(unsorted)
        expected = ["##MuTect", "##Strelka", "##VarScan", "##INFO", "##FORMAT"]
        self.assertEquals(expected, actual)

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


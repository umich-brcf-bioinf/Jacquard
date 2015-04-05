# pylint: disable=line-too-long,too-many-public-methods,invalid-name,too-few-public-methods
from __future__ import print_function, absolute_import, division

import jacquard.variant_callers.zscore_caller as zscore_caller
import jacquard.vcf as vcf
import numpy as np
import test.test_case as test_case
from test.vcf_test import MockVcfReader


class MockTag(object):
    def __init__(self, field_name, sample_values, metaheaders=None):
        self.field_name = field_name
        self.sample_values = sample_values
        if metaheaders:
            self.metaheaders = metaheaders
        else:
            self.metaheaders = []

    def add_tag_values(self, vcf_record):
        vcf_record.add_sample_tag_value(self.field_name, self.sample_values)


def mean(array):
    #pylint: disable=no-member
    return np.mean(array)

def stdev(array):
    #pylint: disable=no-member
    return np.std(array)

class ZScoreTagTest(test_case.JacquardBaseTestCase):
    def test_init_metaheaders(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"13"}, "SB":{"X":"16"}})
        reader = MockVcfReader(records=[rec1, rec2])

        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        self.assertEquals(3, len(tag.metaheaders))
        it = iter(tag.metaheaders)
        self.assertEquals(next(it), '##jacquard.summarize.ZScoreX.X_mean=' + str(tag._mean))
        self.assertEquals(next(it), '##jacquard.summarize.ZScoreX.X_stdev=' + str(tag._stdev))
        self.assertRegexpMatches(next(it), '##FORMAT=<ID=ZScoreX,Number=1,Type=Float,Description="ZScore for X">')

    def test_init_setsPopulationStatistics(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"13"}, "SB":{"X":"16"}})
        reader = MockVcfReader(records=[rec1, rec2])

        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        values = [4, 7, 13, 16]
        self.assertEquals(mean(values), tag._mean)
        self.assertEquals(stdev(values), tag._stdev)


    def test_init_setsPopulationStatisticsParsesFloats(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"2"}, "SB":{"X":"3.5"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"6.5"}, "SB":{"X":"8"}})
        reader = MockVcfReader(records=[rec1, rec2])

        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        values = [2, 3.5, 6.5, 8]
        self.assertEquals(mean(values), tag._mean)
        self.assertEquals(stdev(values), tag._stdev)


    def test_init_setsPopulationStatisticsUsingMaxRangeForMultiValuedInput(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"13,0"}, "SB":{"X":"0,16"}})
        reader = MockVcfReader(records=[rec1, rec2])

        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        values = [4, 7, 13, 16]
        self.assertEquals(mean(values), tag._mean)
        self.assertEquals(stdev(values), tag._stdev)

    def test_init_setsPopulationStatisticsSkipsBlanks(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"."}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"7"}, "SB":{"X":"13"}})
        rec3 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"16"}, "SB":{"X":"."}})
        reader = MockVcfReader(records=[rec1, rec2, rec3])

        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        values = [4, 7, 13, 16]
        self.assertEquals(mean(values), tag._mean)
        self.assertEquals(stdev(values), tag._stdev)

    def test_init_setsPopulationStatisticsConsidersZeros(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"0"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"7"}, "SB":{"X":"13"}})
        rec3 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"16"}, "SB":{"X":"0"}})
        reader = MockVcfReader(records=[rec1, rec2, rec3])

        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        values = [4, 0, 7, 13, 16, 0]
        self.assertAlmostEquals(mean(values), float(tag._mean))
        self.assertAlmostEquals(stdev(values), float(tag._stdev))


    def test_init_setsPopulationStatisticsSkipsSamplesLackingSourceTag(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"Y":"1"}, "SB":{"Y":"2"}})
        rec3 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"13"}, "SB":{"X":"16"}})
        reader = MockVcfReader(records=[rec1, rec2, rec3])

        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        values = [4, 7, 13, 16]
        self.assertEquals(mean(values), tag._mean)
        self.assertEquals(stdev(values), tag._stdev)

    def test_init_setsPopulationStatisticsIgnoresUnparsableValues(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"A1"}, "SB":{"X":"2A"}})
        rec3 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"13"}, "SB":{"X":"16"}})
        reader = MockVcfReader(records=[rec1, rec2, rec3])

        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        values = [4, 7, 13, 16]
        self.assertEquals(mean(values), tag._mean)
        self.assertEquals(stdev(values), tag._stdev)

    def test_init_setsPopulationStatisticsAssignsStddevCorrectlyWhenOneValue(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"."}})
        reader = MockVcfReader(records=[rec1])

        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        values = [4]
        self.assertEquals(mean(values), tag._mean)
        self.assertEquals(stdev(values), tag._stdev)

    def test_init_setsPopulationStatisticsAssignsStddevCorrectlyWhenNoValues(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"."}, "SB":{"X":"."}})
        reader = MockVcfReader(records=[rec1])

        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        expected_mean = None
        expected_stdev = None
        self.assertEquals(expected_mean, tag._mean)
        self.assertEquals(expected_stdev, tag._stdev)


    def test_add_tag(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"13"}, "SB":{"X":"16"}})
        reader = MockVcfReader(records=[rec1, rec2])
        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        tag.add_tag_values(rec1)
        sampleA_tag_values = rec1.sample_tag_values["SA"]
        self.assertEquals("-1.26", sampleA_tag_values["ZScoreX"])
        sampleB_tag_values = rec1.sample_tag_values["SB"]
        self.assertEquals("-0.63", sampleB_tag_values["ZScoreX"])

    def test_add_tag_nullInputsProduceNullZScores(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"."}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"."}, "SB":{"X":"8"}})
        reader = MockVcfReader(records=[rec1, rec2])
        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        tag.add_tag_values(rec1)
        sampleA_tag_values = rec1.sample_tag_values["SA"]
        self.assertEquals("-1.0", sampleA_tag_values["ZScoreX"])
        sampleB_tag_values = rec1.sample_tag_values["SB"]
        self.assertEquals(".", sampleB_tag_values["ZScoreX"])

    def test_add_tag_zeroInputsIncluded(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"0"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"0"}, "SB":{"X":"8"}})
        reader = MockVcfReader(records=[rec1, rec2])
        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        tag.add_tag_values(rec1)
        sampleA_tag_values = rec1.sample_tag_values["SA"]
        self.assertEquals("0.30", sampleA_tag_values["ZScoreX"])
        sampleB_tag_values = rec1.sample_tag_values["SB"]
        self.assertEquals("-0.90", sampleB_tag_values["ZScoreX"])

    def test_add_tag_doesNothingIfNoStdev(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"4"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"}, "SB":{"X":"4"}})
        reader = MockVcfReader(records=[rec1, rec2])
        tag = zscore_caller._ZScoreTag("ZScoreX", "ZScore for X", "X", reader)

        tag.add_tag_values(rec1)

        self.assertEquals(0, tag._stdev)
        self.assertEqual(["X"], sorted(rec1.sample_tag_values["SA"].keys()))
        self.assertEqual(["X"], sorted(rec1.sample_tag_values["SB"].keys()))

class AlleleFreqZScoreTagTest(test_case.JacquardBaseTestCase):
    def test_init_metaheaders(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_AF_RANGE":"4"},
                                                "SB":{"JQ_SUMMARY_AF_RANGE":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_AF_RANGE":"13"},
                                                "SB":{"JQ_SUMMARY_AF_RANGE":"16"}})
        reader = MockVcfReader(records=[rec1, rec2])

        tag = zscore_caller._AlleleFreqZScoreTag(reader)

        self.assertEquals(3, len(tag.metaheaders))
        it = iter(tag.metaheaders)
        self.assertRegexpMatches(next(it), '##jacquard.summarize.JQ_SUMMARY_AF_ZSCORE.JQ_SUMMARY_AF_RANGE_mean=')
        self.assertRegexpMatches(next(it), '##jacquard.summarize.JQ_SUMMARY_AF_ZSCORE.JQ_SUMMARY_AF_RANGE_stdev=')
        self.assertRegexpMatches(next(it), '##FORMAT=<ID=JQ_SUMMARY_AF_ZSCORE,Number=1,Type=Float,Description="Concordance of reported allele frequencies.*">')

    def test_add_tag(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_AF_RANGE":"4"},
                                                "SB":{"JQ_SUMMARY_AF_RANGE":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_AF_RANGE":"13"},
                                                "SB":{"JQ_SUMMARY_AF_RANGE":"16"}})
        reader = MockVcfReader(records=[rec1, rec2])
        tag = zscore_caller._AlleleFreqZScoreTag(reader)

        tag.add_tag_values(rec1)
        sampleA_tag_values = rec1.sample_tag_values["SA"]
        self.assertEquals("-1.26", sampleA_tag_values["JQ_SUMMARY_AF_ZSCORE"])
        sampleB_tag_values = rec1.sample_tag_values["SB"]
        self.assertEquals("-0.63", sampleB_tag_values["JQ_SUMMARY_AF_ZSCORE"])

class DepthZScoreTagTest(test_case.JacquardBaseTestCase):
    def test_init_metaheaders(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_DP_RANGE":"4"},
                                                "SB":{"JQ_SUMMARY_DP_RANGE":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_DP_RANGE":"13"},
                                                "SB":{"JQ_SUMMARY_DP_RANGE":"16"}})
        reader = MockVcfReader(records=[rec1, rec2])

        tag = zscore_caller._DepthZScoreTag(reader)

        self.assertEquals(3, len(tag.metaheaders))
        it = iter(tag.metaheaders)
        self.assertRegexpMatches(next(it), '##jacquard.summarize.JQ_SUMMARY_DP_ZSCORE.JQ_SUMMARY_DP_RANGE_mean=')
        self.assertRegexpMatches(next(it), '##jacquard.summarize.JQ_SUMMARY_DP_ZSCORE.JQ_SUMMARY_DP_RANGE_stdev=')
        self.assertRegexpMatches(next(it), '##FORMAT=<ID=JQ_SUMMARY_DP_ZSCORE,Number=1,Type=Float,Description="Concordance of reported depth.*">')

    def test_add_tag(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_DP_RANGE":"4"},
                                                "SB":{"JQ_SUMMARY_DP_RANGE":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_DP_RANGE":"13"},
                                                "SB":{"JQ_SUMMARY_DP_RANGE":"16"}})
        reader = MockVcfReader(records=[rec1, rec2])
        tag = zscore_caller._DepthZScoreTag(reader)

        tag.add_tag_values(rec1)
        sampleA_tag_values = rec1.sample_tag_values["SA"]
        self.assertEquals("-1.26", sampleA_tag_values["JQ_SUMMARY_DP_ZSCORE"])
        sampleB_tag_values = rec1.sample_tag_values["SB"]
        self.assertEquals("-0.63", sampleB_tag_values["JQ_SUMMARY_DP_ZSCORE"])


class ZScoreCallerTest(test_case.JacquardBaseTestCase):
    def test_init_createsAllTags(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_AF_RANGE":"4"},
                                                "SB":{"JQ_SUMMARY_AF_RANGE":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_AF_RANGE":"13"},
                                                "SB":{"JQ_SUMMARY_AF_RANGE":"16"}})
        reader = MockVcfReader(records=[rec1, rec2])

        caller = zscore_caller.ZScoreCaller(reader)

        self.assertEquals(2, len(caller._tags))
        tag_ids = [tag.tag_id for tag in caller._tags]
        self.assertIn("JQ_SUMMARY_AF_ZSCORE", tag_ids)
        self.assertIn("JQ_SUMMARY_DP_ZSCORE", tag_ids)

    def test_init_createsAllMetaheaders(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_AF_RANGE":"4"},
                                                "SB":{"JQ_SUMMARY_AF_RANGE":"7"}})
        rec2 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"JQ_SUMMARY_AF_RANGE":"13"},
                                                "SB":{"JQ_SUMMARY_AF_RANGE":"16"}})
        reader = MockVcfReader(records=[rec1, rec2])

        caller = zscore_caller.ZScoreCaller(reader)

        self.assertEquals(6, len(caller.metaheaders))
        it = iter(caller.metaheaders)
        self.assertRegexpMatches(next(it), "JQ_SUMMARY_AF_ZSCORE.JQ_SUMMARY_AF_RANGE_mean")
        self.assertRegexpMatches(next(it), "JQ_SUMMARY_AF_ZSCORE.JQ_SUMMARY_AF_RANGE_stdev")
        self.assertRegexpMatches(next(it), "FORMAT=<ID=JQ_SUMMARY_AF_ZSCORE")
        self.assertRegexpMatches(next(it), "JQ_SUMMARY_DP_ZSCORE.JQ_SUMMARY_DP_RANGE_mean")
        self.assertRegexpMatches(next(it), "JQ_SUMMARY_DP_ZSCORE.JQ_SUMMARY_DP_RANGE_stdev")
        self.assertRegexpMatches(next(it), "FORMAT=<ID=JQ_SUMMARY_DP_ZSCORE")

    def test_add_tags(self):
        rec1 = vcf.VcfRecord("1", "42", "A", "C",
                             sample_tag_values={"SA":{"X":"4"},
                                                "SB":{"X":"7"}})
        reader = MockVcfReader(records=[rec1])

        caller = zscore_caller.ZScoreCaller(reader)
        caller._tags = [MockTag("Y", {"SA":"A42", "SB":"B42"})]
        caller.add_tags(rec1)

        self.assertEquals({"X":"4", "Y":"A42"}, rec1.sample_tag_values["SA"])

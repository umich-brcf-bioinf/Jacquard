#!/usr/bin/python2.7
import pandas as pd
from pandas import *
import unittest
from pandas.util.testing import assert_frame_equal
import pandas.util.testing as tm
from bin.rollup import gene_rollup_highest_impact, gene_rollup_damaging_impact, combine_dfs
from StringIO import StringIO
import pprint
import os
from os import listdir
from os.path import isfile, join

pd.set_option('chained_assignment', None)

def dataframe(input_data, sep="\t", index_col=None):
    return pd.read_csv(StringIO(input_data), sep=sep, header=False, dtype='str', index_col=index_col)
    
class RollupTestCase(unittest.TestCase):
    def test_gene_rollup_highest_impact(self):
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	HIGHEST_IMPACT	RNA_FREQ_sampleA	RNA_FREQ_sampleB
sample1	1	2	A	T	foo	HIGH	0.5	0.2
sample1	1	2	A	G	foo	LOW	0.7	.28
sample2	2	3	A	T	bar	MODIFIER	0.5	0.1
sample2	2	3	A	T	foo	.	0.5	0.1
sample2	2	3	A	T	bar	MODERATE	0.5	.'''
        df = dataframe(input_string)
        
        actual_df = gene_rollup_highest_impact(df)
 
        impact_RNA_FREQ_sampleA = pd.Series({"bar": "mx", "foo": "hl"}, name="GENE_SYMBOL", index=["bar", "foo"])
        impact_RNA_FREQ_sampleB = pd.Series({"bar": "x", "foo": "hl"}, name="GENE_SYMBOL", index=["bar", "foo"])
        
        score_RNA_FREQ_sampleA = pd.Series({"bar": 1.00000, "foo": 100000.00001}, name="GENE_SYMBOL", index=["bar", "foo"])
        score_RNA_FREQ_sampleB = pd.Series({"bar": 0.00000, "foo": 100000.00001}, name="GENE_SYMBOL", index=["bar", "foo"])

        tm.assert_series_equal(impact_RNA_FREQ_sampleA, actual_df["Impact_complete", "RNA_FREQ_sampleA"])
        tm.assert_series_equal(impact_RNA_FREQ_sampleB, actual_df["Impact_complete", "RNA_FREQ_sampleB"])
       
        tm.assert_series_equal(score_RNA_FREQ_sampleA, actual_df["Impact_score", "RNA_FREQ_sampleA"])
        tm.assert_series_equal(score_RNA_FREQ_sampleB, actual_df["Impact_score", "RNA_FREQ_sampleB"])
        
    def test_gene_rollup_damaging_impact(self): 
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	Impact_Damaging	RNA_FREQ_sampleA	RNA_FREQ_sampleB
sample1	1	2	A	T	foo	0	0.5	0.2
sample1	1	2	A	G	foo	1	0.7	.28
sample2	2	3	A	T	bar	0	.	0.1
sample2	2	3	A	T	foo	7	0.5	0.1
sample2	2	3	A	T	bar	2	0.5	.'''
        df = dataframe(input_string)
        
        actual_df = gene_rollup_damaging_impact(df)
        
        impact_RNA_FREQ_sampleA = pd.Series({"bar": 2, "foo": 8}, name="GENE_SYMBOL", index=["bar", "foo"])
        impact_RNA_FREQ_sampleB = pd.Series({"bar": 2, "foo": 8}, name="GENE_SYMBOL", index=["bar", "foo"])
        
        tm.assert_series_equal(impact_RNA_FREQ_sampleA, actual_df["Impact_Damaging", "RNA_FREQ_sampleA"])
        tm.assert_series_equal(impact_RNA_FREQ_sampleB, actual_df["Impact_Damaging", "RNA_FREQ_sampleB"])
        
        
    def test_combine_dfs(self):
        input_string = \
'''SAMPLE_NAME	CHROM	POS	REF	ANNOTATED_ALLELE	GENE_SYMBOL	HIGHEST_IMPACT	Impact_Damaging	RNA_FREQ_sampleA	RNA_FREQ_sampleB
sample1	1	2	A	T	foo	HIGH	0	0.5	0.2
sample1	1	2	A	G	foo	MODERATE	1	0.7	.28
sample2	2	3	A	T	bar	LOW	0	.	0.1
sample2	2	3	A	T	foo	MODIFIER	7	0.5	0.1
sample2	2	3	A	T	bar	LOW	2	0.5	.'''
        df1 = dataframe(input_string)
        df2 = dataframe(input_string)
        
        df1 = gene_rollup_highest_impact(df1)
        df2 = gene_rollup_damaging_impact(df2)
        
        combined_df = combine_dfs(df1, df2)

        impact_combined_RNA_FREQ_sampleA = pd.Series({"bar": "l", "foo": "hmx"}, name="GENE_SYMBOL", index=["foo", "bar"])
        impact_combined_RNA_FREQ_sampleB = pd.Series({"bar": "l", "foo": "hmx"}, name="GENE_SYMBOL", index=["foo", "bar"])
        
        impact_damaging_RNA_FREQ_sampleA = pd.Series({"bar": 2, "foo": 8}, name="GENE_SYMBOL", index=["foo", "bar"])
        impact_damaging_RNA_FREQ_sampleB = pd.Series({"bar": 2, "foo": 8}, name="GENE_SYMBOL", index=["foo", "bar"])
        
        rank = pd.Series({"bar": 2., "foo": 1.}, name="GENE_SYMBOL", index=["foo", "bar"])
        
        tm.assert_series_equal(impact_combined_RNA_FREQ_sampleA, combined_df["Impact_complete", "RNA_FREQ_sampleA"])
        tm.assert_series_equal(impact_combined_RNA_FREQ_sampleB, combined_df["Impact_complete", "RNA_FREQ_sampleB"])

        tm.assert_series_equal(impact_damaging_RNA_FREQ_sampleA, combined_df["Impact_Damaging", "RNA_FREQ_sampleA"])
        tm.assert_series_equal(impact_damaging_RNA_FREQ_sampleB, combined_df["Impact_Damaging", "RNA_FREQ_sampleB"])
        
        tm.assert_series_equal(rank, combined_df["Impact_complete_Rank"])
        tm.assert_series_equal(rank, combined_df["Impact_Damaging_Rank"])
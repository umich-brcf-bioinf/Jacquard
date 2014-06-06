#!/usr/bin/python2.7
import math
import numpy
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
from pandas import *
import random
import re
import sys 
from StringIO import StringIO


dataString1 = \
'''CHROM	POS	REF	ALT	SAMPLE_NAME	GT	DP
chr1	1	A	T	sample1	1/1	12.8%
chr1	2	A	T	sample1	0/1	2.4
chr1	3	A	T	sample2	1/1	0.03
chr1	4	A	T	sample3	0/1	4'''

df = pd.read_csv(StringIO(dataString1), sep="\t", header=False, dtype='str')

print df.applymap(lambda x: type(x))
# writer = ExcelWriter("test_to_excel.xlsx")
# df.to_excel(writer, "Variant_output", index=True, merge_cells = 0)  
# writer.save() 

df.to_csv("test_to_excel.csv", index=True)
    
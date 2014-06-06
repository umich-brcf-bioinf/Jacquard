#!/usr/bin/python2.7
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
from pandas import *
import sys 
from StringIO import StringIO

def dataframe(input_data, sep="\t", index_col=None):
    return pd.read_csv(StringIO(input_data), sep=sep, header=False, dtype='str', index_col=index_col)
    
script_dir = os.path.dirname(os.path.abspath(__file__))

expected_string = \
'''\tSAMPLE_NAME\tFORMAT\tSample1\tGT
0	sample1	GT:DP	A:1	A
1	sample2	GT:DP	B:2	B
2	sample3	DP	3	NaN
3	sample6	GT	D	D'''
df = dataframe(expected_string, index_col=0)  
   
writer = ExcelWriter(script_dir + '/test_output.xlsx')
    
df.to_excel(writer, "Variant_output", index=True, merge_cells = 0)  
writer.save()    
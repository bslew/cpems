#!/usr/bin/env python

import sys
#from pyCPEDScommonFunctions import cpedsPythCommon

import numpy as np
import os
import sys
from optparse import OptionParser


programDescription="""
Calculates mean function value from input function files 

Example usage:
average_fn.py file1 file2 ...


fileX - input file with 2 columns

"""

parser = OptionParser(description=programDescription)
# parser.add_option("", "--test", action="store_true", dest="test", default=False, help="DEBUG.")
parser.add_option("", "--stdev", action="store_true", dest="stdev", default=False, help="calculate also stdev")
parser.add_option("-o", "--output", dest="output", default='', type="string", help='output file name', metavar="VALUE")
parser.add_option("-x", "--colX", dest="colX", default=0, type="int", help='column to be used for averaging (0)', metavar="VALUE")
parser.add_option("-y", "--colY", dest="colY", default=1, type="int", help='column to be used for averaging (1)', metavar="VALUE")
# parser.add_option("", "--sortX", action="store_true", dest="sortX", help='sort output by X values')

(option, args) = parser.parse_args()

sum=0
N=len(args)
stdev=0

for i in range(N):
    d=np.loadtxt(args[i])
    d=d[:,[option.colX,option.colY]]
    if i==0:
        sum=d
        if option.stdev:
            stdev=d[:,1]
    else:
        sum=sum+d
        if option.stdev:
            stdev=np.vstack([stdev,d[:,1]])
        
out=sum/N
if option.stdev:
    stdev=np.std(stdev, 0)
    stdev=np.reshape(stdev, [len(stdev),1])
#     print stdev
#     print sum
    out=np.hstack([sum/N,stdev])

# out=out.sort(key=lambda x: x[1])
if option.output!='':
    np.savetxt(option.output,out)
else:
#     print out
    np.savetxt(sys.stdout,out)

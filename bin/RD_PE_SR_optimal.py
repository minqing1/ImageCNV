#! /usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import cnvnator_filter
import lumpy_filter
import merge_cnvsv
import cnv_filter 
import pandas as pd
from datetime import datetime
class Run(object):
    def __init__(self):
        self.tmp = '../tmp/'
    def __Define__(self,samplename):
        self.samplename = samplename

    def main(self):
        usage = 'python run.py wgslist.txt,error,please check!'
        if len(sys.argv)<2:
            print usage
            sys.exit()

        else:
            # single wgs sample
            samplename = sys.argv[1]
            cnv = cnvnator_filter.cnvnatorFilter(samplename)
            cnv.main()
            sv = lumpy_filter.lumpyFilter(samplename)
            sv.main()
            merge = merge_cnvsv.cnvMerged(samplename)
            merge.merged()
            cnv_filter.hundredFilter(samplename)
if __name__ == '__main__':
    run = Run()
    run.main()

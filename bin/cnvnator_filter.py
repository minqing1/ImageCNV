#! /usr/bin/env python
# -*- coding:utf-8 -*-
# Desciption:RD cnv set filter hg19 gap and chr6 MHC

import os
import sys
from datetime import datetime

class cnvnatorFilter(object):
    def __init__(self,samplename):
        self.chrom_list = ['chr'+str(x) for x in range(1,23)]
        self.samplename = samplename
        self.gap_path ='database/gap/'
        self.segdup_path = 'database/seg_dup/'
        self.precnv_path = 'example/'
        self.hcbed = 'database/ceph18_high_coverage/'
    def qFilter(self):
        '''filter q0>0.5 and retain length >= 1kb cnv'''
        for chrom in self.chrom_list:
            os.system('grep -w '+chrom+' '+self.precnv_path+self.samplename+'.cnv.xls > '+self.precnv_path+chrom+'.cnv.xls')
            fileopen=open(self.precnv_path+chrom+'.cnv.xls','r')
            filesave = open(self.precnv_path+chrom+'.cnv.qFilter.xls','w')
            for i in fileopen:
                cnvlist = i.strip('\n').split()
                q_num = float(cnvlist[-1])
                info = cnvlist[3].split(':')[-1]
                start = int(info.split('-')[0])
                end = int(info.split('-')[1])
                length = end - start
                if q_num <= 0.5 and length>=1000:
                #if length >= 1000:
                    filesave.write(i)
            fileopen.close()
            filesave.close()
        print 'q0 is done'
             
    def gapFilter(self):
        ''' filter hg19 gap '''
        for chrom in self.chrom_list:
            fileopen = open(self.precnv_path+chrom+'.cnv.qFilter.xls','r')
            filesave2 = open(self.precnv_path+chrom+'.cnv.in_gap.xls','w')
            # hg19gap overlap>50%
            for i in fileopen:
                cnvlist = i.split()
                info = cnvlist[3]
                postion = info.split(':')[1]
                cnv_start = int(postion.split('-')[0])
                cnv_end = int(postion.split('-')[1])
                gap = open(self.gap_path+chrom+'_gap.txt','r')
                for j in gap:
                    gaplist = j.split()
                    gap_start = int(gaplist[2])
                    gap_end = int(gaplist[3])
                    #1.back overlap in cnv >50%
                    if cnv_start <= gap_start< cnv_end <= gap_end:
                        overlap_length = cnv_end - gap_start
                        cnv_length = cnv_end - cnv_start
                        overlap = overlap_length / float(cnv_length)
                        if overlap >= 0.7:
                            filesave2.write(i)
                    #2.forehead overlap in cnv >50%
                    elif gap_start <= cnv_start < gap_end <= cnv_end:
                        overlap_length = gap_end - cnv_start
                        cnv_length = cnv_end - cnv_start
                        overlap = overlap_length / float(cnv_length)
                        if overlap >= 0.7:
                            filesave2.write(i)
                    #3. gap include in cnv 100%
                    elif cnv_start < gap_start < gap_end < cnv_end:
                        overlap_length = gap_end - gap_start
                        cnv_length = cnv_end - cnv_start
                        overlap = overlap_length / float(cnv_length)
                        if overlap >= 0.7:
                            filesave2.write(i)
                    #4 cnv == gap
                    elif cnv_start == gap_start and cnv_end == gap_end:
                        filesave2.write(i)
                    elif gap_start < cnv_start < cnv_end < gap_end:
                        filesave2.write(i)
                    else:
                        continue
                gap.close()
            #different row in filtered.xls and cnv.xls
            fileopen.close()
            filesave2.close()
            os.system('grep -xvf '+self.precnv_path+chrom+'.cnv.in_gap.xls '+self.precnv_path+chrom+'.cnv.qFilter.xls >'+self.precnv_path+chrom+'.cnv.filter_gap.xls')
        print 'filter hg19 gap is done'
    def segdupFilter(self):
        ''' filter hg19 segdup 
            inputfile : NA12878.cnv.chrom.filter_gap.xls
        '''
        for chrom in self.chrom_list:
            fileopen = open(self.precnv_path+chrom+'.cnv.filter_gap.xls','r')
            filesave2 = open(self.precnv_path+chrom+'.cnv.in_segdup.xls','w')
            for i in fileopen:
                cnvlist = i.split()
                #q0 argument
                #Q_num = float(cnvlist[-1])
                #if Q_num > 0.5:
                    #continue
                info = cnvlist[3]
                postion = info.split(':')[1]
                cnv_start = int(postion.split('-')[0])
                cnv_end = int(postion.split('-')[1])
                dup = open(self.segdup_path+chrom+'_dup.txt','r')
                for j in dup:
                    duplist = j.split()
                    dup_chr = duplist[0]
                    dup_start = int(duplist[1])
                    dup_end = int(duplist[2])
                    if cnv_start <= dup_start < cnv_end <= dup_end:
                        overlap_length = cnv_end - dup_start
                        cnv_length = cnv_end - cnv_start
                        overlap = overlap_length / float(cnv_length)
                        if overlap >= 0.7:
                            filesave2.write(i)
                    #2.forehead overlap in cnv >50%
                    elif dup_start <= cnv_start < dup_end <= cnv_end:
                        overlap_length = dup_end - cnv_start
                        cnv_length = cnv_end - cnv_start
                        overlap = overlap_length / float(cnv_length)
                        if overlap >= 0.7:
                            filesave2.write(i)
                    #3. dup include in cnv 100%
                    elif cnv_start < dup_start < dup_end < cnv_end:
                        overlap_length = dup_end - dup_start
                        cnv_length = cnv_end - cnv_start
                        overlap = overlap_length / float(cnv_length)
                        if overlap >= 0.7:
                            filesave2.write(i)
                    #4 cnv == dup
                    elif cnv_start == dup_start and cnv_end == dup_end:
                        filesave2.write(i)
                    #5 cnv include in dup 100%
                    elif dup_start < cnv_start < cnv_end < dup_end:
                        filesave2.write(i)
                    elif cnv_start < dup_start < dup_end < cnv_end:
                        filesave2.write(i)
                    else:
                        continue
                dup.close()
            fileopen.close()
            filesave2.close()
            os.system('grep -xvf '+self.precnv_path+chrom+'.cnv.in_segdup.xls '+self.precnv_path+chrom+'.cnv.filter_gap.xls >'+self.precnv_path+chrom+'.cnv.filter_dup.xls')
        print 'filter hg19 segdup is done'
    def filterHc(self):
        """filter high coverage region"""
        for chrom in self.chrom_list:
            fileopen = open(self.precnv_path+chrom+'.cnv.filter_dup.xls','r')
            filesave2 = open(self.precnv_path+chrom+'.cnv.in_hc.xls','w')
            for i in fileopen:
                cnvlist = i.split()
                #q0 argument
                #Q_num = float(cnvlist[-1])
                #if Q_num > 0.5:
                    #continue
                info = cnvlist[3]
                postion = info.split(':')[1]
                cnv_start = int(postion.split('-')[0])
                cnv_end = int(postion.split('-')[1])
                dup = open(self.hcbed+chrom+'_hc.txt','r')
                for j in dup:
                    duplist = j.split()
                    dup_chr = duplist[0]
                    dup_start = int(duplist[1])
                    dup_end = int(duplist[2])
                    if cnv_start <= dup_start < cnv_end <= dup_end:
                        overlap_length = cnv_end - dup_start
                        cnv_length = cnv_end - cnv_start
                        overlap = overlap_length / float(cnv_length)
                        if overlap >= 0.1:
                            filesave2.write(i)
                    #2.forehead overlap in cnv >50%
                    elif dup_start <= cnv_start < dup_end <= cnv_end:
                        overlap_length = dup_end - cnv_start
                        cnv_length = cnv_end - cnv_start
                        overlap = overlap_length / float(cnv_length)
                        if overlap >= 0.1:
                            filesave2.write(i)
                    #3. dup include in cnv 100%
                    elif cnv_start < dup_start < dup_end < cnv_end:
                        overlap_length = dup_end - dup_start
                        cnv_length = cnv_end - cnv_start
                        overlap = overlap_length / float(cnv_length)
                        if overlap >= 0.1:
                            filesave2.write(i)
                    #4 cnv == dup
                    elif cnv_start == dup_start and cnv_end == dup_end:
                        filesave2.write(i)
                    #5 cnv include in dup 100%
                    elif dup_start < cnv_start < cnv_end < dup_end:
                        filesave2.write(i)
                    elif cnv_start < dup_start < dup_end < cnv_end:
                        filesave2.write(i)
                    else:
                        continue
                dup.close()
            fileopen.close()
            filesave2.close()
            os.system('grep -xvf '+self.precnv_path+chrom+'.cnv.in_hc.xls '+self.precnv_path+chrom+'.cnv.filter_dup.xls >'+self.precnv_path+chrom+'.cnv.final.xls')
        os.system('cat '+self.precnv_path+'chr*.cnv.final.xls > '+self.precnv_path+self.samplename+'.cnv.final.xls')
        print 'filter high coverage region is done'
    def compare_gscnv(self):
        for chrom in self.chrom_list:
            #fileopen = open(self.precnv_path+chrom+'.cnv.final.xls','r')
            fileopen = open(self.precnv_path+chrom+'.cnv.qFilter.xls','r')
            filesave1 = open(self.precnv_path+chrom+'.cnv.50.tmp.xls','w')
            filesave2 = open(self.precnv_path+chrom+'.cnv.1050.tmp.xls','w')
            # gscnv overlap
            for i in fileopen:
                cnvlist = i.split()
                info = cnvlist[3]
                postion = info.split(':')[1]
                cnv_start = int(postion.split('-')[0])
                cnv_end = int(postion.split('-')[1])
                cnv_type = cnvlist[2]
                gscnv = open(self.gscnv_path+chrom+'_GSCNV_1000.txt','r')
                #gscnv = open(self.gsdeletion+chrom+'_deletion.txt','r')
                #gscnv = open(self.Mills+chrom+'_deletion.txt','r')
                for j in gscnv:
                    gaplist = j.split()
                    gap_start = int(gaplist[1])
                    gap_end = int(gaplist[2])
                    gap_type = gaplist[4].lower()
                    if gap_type in cnv_type:
                    #if cnv_type == 'deletion':
                    #1.back overlap in cnv >50%
                        if cnv_start <= gap_start and gap_start< cnv_end <= gap_end:
                            overlap_length = cnv_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                filesave1.write(i)
                            elif 0.1<=overlap_cnv<0.5 and 0.1<= overlap_gap<0.5:
                                filesave2.write(i)
                            elif 0.1<=overlap_cnv<0.5 and 0.5<= overlap_gap:
                                filesave2.write(i)
                            elif 0.5<=overlap_cnv and 0.1<= overlap_gap<0.5:
                                filesave2.write(i)
                        #2.forehead overlap in cnv >50%
                        elif cnv_end >= gap_end and gap_start <= cnv_start < gap_end:
                            overlap_length = gap_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                filesave1.write(i)
                            elif 0.1 <= overlap_cnv<0.5 and 0.1 <= overlap_gap <0.5:
                                filesave2.write(i)
                            elif 0.1<=overlap_cnv<0.5 and 0.5<= overlap_gap:
                                filesave2.write(i)
                            elif 0.5<=overlap_cnv and 0.1<= overlap_gap<0.5:
                                filesave2.write(i)
                        #3. gap include in cnv 100%
                        elif cnv_start < gap_start < gap_end < cnv_end:
                            overlap_length = gap_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                filesave1.write(i)
                            elif 0.1<=overlap_cnv<0.5 and 0.1<=overlap_gap<0.5:
                                filesave2.write(i)
                            elif 0.1<=overlap_cnv<0.5 and 0.5<= overlap_gap:
                                filesave2.write(i)
                            elif 0.5<=overlap_cnv and 0.1<= overlap_gap<0.5:
                                filesave2.write(i)
                        #4 cnv == gap
                        elif cnv_start == gap_start and cnv_end == gap_end:
                            filesave1.write(i)
                        #5 cnv in gap
                        elif gap_start < cnv_start < cnv_end <gap_end:
                            overlap_length = cnv_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                filesave1.write(i)
                            elif 0.1<=overlap_cnv<0.5 and 0.1<=overlap_gap<0.5:
                                filesave2.write(i)
                            elif 0.1<=overlap_cnv<0.5 and 0.5<= overlap_gap:
                                filesave2.write(i)
                            elif 0.5<=overlap_cnv and 0.1<= overlap_gap<0.5:
                                filesave2.write(i)
                        else:
                            continue
                gscnv.close()
            #different row in filtered.xls and cnv.xls
            fileopen.close()
            filesave1.close()
            filesave2.close()
            os.system('sort -u '+self.precnv_path+chrom+'.cnv.50.tmp.xls > '+self.precnv_path+chrom+'.cnv.50.xls')
            os.system('sort -u '+self.precnv_path+chrom+'.cnv.1050.tmp.xls > '+self.precnv_path+chrom+'.cnv.1050.xls')
            
    def main(self):
        self.qFilter()
        self.gapFilter()
        self.segdupFilter()
        self.filterHc()
        #run2 = self.compare_gscnv()
        print 'cnvnator filter is done at %s' % datetime.now()
if __name__ == '__main__':
    print 'start at %s' % datetime.now()
    import sys
    usage = 'python cnvnator_filter.py binsize coverage'
    if len(sys.argv)<2:
        print usage
        sys.exit()
    binsize = str(sys.argv[1])
    coverage = str(sys.argv[2])
    Run = cnvnatorFilter(binsize,coverage)
    run = Run.main()

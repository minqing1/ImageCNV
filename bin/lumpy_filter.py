#! /usr/bin/env python
# -*- coding:utf-8 -*-
# Description:lumpy sv filter script 
# 1.filter hg 19 gap  2.filter hg19 segdup 3.filter readcount
import os
import sys
from datetime import datetime

class lumpyFilter(object):
    def __init__(self,samplename):
        self.samplename = samplename
        self.path = 'example/'
        self.out = 'example/'+samplename
        self.chromlist = ['chr'+str(x) for x in range(1,23)]
        self.gap_path = 'database/gap/'
        self.segdup_path = 'database/seg_dup/'
        self.hcbed = 'database/ceph18_high_coverage/'
    def filterType(self):
        """ only retain dup and del
        """
        f = open(self.path+self.samplename+'.sv.result.xls','r')
        fsave = open(self.out+'.sv.dupdel.result.xls','w')
        for i in f:
            linelist = i.split()
            Type = linelist[14]
            # retain dup and del
            if Type == 'TYPE:INTERCHROM' or Type == 'TYPE:INVERSION' or Type == 'type':
                continue
            else:
                fsave.write(i)
        fsave.close()
        print 'retain dup and del is done'
    def filterCount(self):
        """ retain length >= 1000 cnv
            inputfile:chrom.sv.dupdel.result.xls
        """
        for chrom in self.chromlist:
            target = self.out+'.'+chrom+'.sv.dupdel.result.xls'
            os.system('grep -w '+chrom+' '+self.out+'.sv.dupdel.result.xls > '+target)
            #c = 24
            f = open(target,'r')
            fsave = open(self.out+'.'+chrom+'.sv.dupdel.tmp.result.xls','w')
            for i in f:
                cnvlist = i.split()
                cnvinfo = cnvlist[17].split(':')
                cnv_start = int(cnvinfo[2].split(';')[0])
                cnv_end = int(cnvinfo[3])
                length = cnv_end - cnv_start
                strand = int(cnvlist[16].split(',')[-1])

                if 150000>=length >= 1000:
                    fsave.write(i)
            f.close()
            fsave.close()
        print 'filter size and read counts is done'
    def sv_classify(self):
        """classify the cnv by evidence,pe,sr,pe+sr"""
        for chrom in self.chromlist:
            #target = self.out+'.'+chrom+'.sv.dupdel.tmp.result.xls'
            target = self.out+'.'+chrom+'.sv.dupdel.final.result.xls'
            f = open(target,'r')
            filesave1 = open(self.out+'.'+chrom+'.sv.dupdel.pe.result.xls','w')
            filesave2 = open(self.out+'.'+chrom+'.sv.dupdel.sr.result.xls','w')
            filesave3 = open(self.out+'.'+chrom+'.sv.dupdel.pe_sr.result.xls','w')
            for i in f:
                cnvlist = i.split()
                evidence = cnvlist[15]
                if evidence.startswith('IDS:1') and ';2' not in evidence:
                    filesave1.write(i)
                elif evidence.startswith('IDS:2'):
                    filesave2.write(i)
                else:
                    filesave3.write(i)
            filesave1.close()
            filesave2.close()
            filesave3.close()
        print 'sv classify is done' 
    def filterRefgap(self):
        """ filter hg19 gap
        """
        for chrom in self.chromlist:
            target = self.out+'.'+chrom+'.sv.dupdel.tmp.result.xls'
            f = open(target,'r')
            filesave2 = open(self.out+'.'+chrom+'.sv.dupdel.in_gap.result.xls','w')
            for i in f:
                cnvlist = i.split()
                cnvinfo = cnvlist[17].split(':')
                cnv_start = int(cnvinfo[2].split(';')[0])
                cnv_end = int(cnvinfo[3])
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
                    elif  gap_start <= cnv_start < gap_end <= cnv_end:
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
            f.close()
            filesave2.close()
            os.system('grep -xvf '+self.out+'.'+chrom+'.sv.dupdel.in_gap.result.xls '+target+' >'+self.out+'.'+chrom+'.sv.dupdel.filter_gap.result.xls')
        #os.system('cat '+self.out+'.chr*'+'.sv.dupdel.pass.result.xls >'+self.out+'.sv.dupdel.pass.result.xls')
        print 'filter hg19 gap is done'
    def filterSegdup(self):
        """ filter hg19 segdup
            inputfile :chrom.sv.dupdel.filter_gap.result.xls
        """
        for chrom in self.chromlist:
            target = self.out+'.'+chrom+'.sv.dupdel.filter_gap.result.xls' 
            fileopen = open(target,'r')
            filesave2 = open(self.out+'.'+chrom+'.sv.dupdel.in_segdup.result.xls','w')
            for i in fileopen:
                cnvlist = i.split()
                cnvinfo = cnvlist[17].split(':')
                cnv_start = int(cnvinfo[2].split(';')[0])
                cnv_end = int(cnvinfo[3])
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
                    else:
                        continue
                dup.close()
            fileopen.close()
            filesave2.close()
            os.system('grep -xvf '+self.out+'.'+chrom+'.sv.dupdel.in_segdup.result.xls '+target+' >'+self.out+'.'+chrom+'.sv.dupdel.filter_dup.result.xls')
        print 'filter hg19 segdup is done'
    def filterHc(self):
        """ filter ceph18 high coverage region """
        for chrom in self.chromlist:
            target = self.out+'.'+chrom+'.sv.dupdel.filter_dup.result.xls' 
            fileopen = open(target,'r')
            filesave2 = open(self.out+'.'+chrom+'.sv.dupdel.in_hc.result.xls','w')
            for i in fileopen:
                cnvlist = i.split()
                cnvinfo = cnvlist[17].split(':')
                cnv_start = int(cnvinfo[2].split(';')[0])
                cnv_end = int(cnvinfo[3])
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
                    else:
                        continue
                dup.close()
            fileopen.close()
            filesave2.close()
            os.system('grep -xvf '+self.out+'.'+chrom+'.sv.dupdel.in_hc.result.xls '+target+' >'+self.out+'.'+chrom+'.sv.dupdel.final.result.xls')
        os.system('cat '+self.out+'.chr*.sv.dupdel.final.result.xls > '+self.out+'.sv.final.result.xls')
        print 'filter hc region is done'
    def compare_gs(self):
        """compare with NA12878 gs cnv set,total 1991 autosome cnv"""
        for chrom in self.chromlist:
            fileopen = open(self.out+'.'+chrom+'.sv.dupdel.final.result.xls','r')
            #fileopen = open(self.out+'.'+chrom+'.sv.dupdel.pe.result.xls','r')
            filesave1 = open(self.out+'.'+chrom+'.sv.merged50.tmp.xls','w')
            filesave2 = open(self.out+'.'+chrom+'.sv.merged1050.tmp.xls','w')
            # gscnv overlap
            for i in fileopen:
                cnvlist = i.split()
                cnvinfo = cnvlist[17].split(':')
                cnv_start = int(cnvinfo[2].split(';')[0])
                cnv_end = int(cnvinfo[3])
                cnv_type = cnvlist[14].split(':')[1]
                gscnv = open(self.gscnv_path+chrom+'_GSCNV_1000.txt','r')
                for j in gscnv:
                    gaplist = j.split()
                    gap_start = int(gaplist[1])
                    gap_end = int(gaplist[2])
                    gap_type = gaplist[4]
                    #1.back overlap in cnv >50%
                    if gap_type in cnv_type:
                        if cnv_start <= gap_start and gap_start< cnv_end <= gap_end:
                            overlap_length = cnv_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end -gap_start
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
                        #2.forehead overlap in cnv >50%
                        elif cnv_end >= gap_end and gap_start <= cnv_start < gap_end:
                            overlap_length = gap_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end -gap_start
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
                        #3. gap include in cnv 100%
                        elif cnv_start < gap_start < gap_end < cnv_end:
                            overlap_length = gap_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end -gap_start
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
                        #4 cnv == gap
                        elif cnv_start == gap_start and cnv_end == gap_end:
                            filesave1.write(i)
                        #5 cnv in gap 100%
                        elif gap_start < cnv_start < cnv_end < gap_end:
                            overlap_length = cnv_end -cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end -gap_start
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
                        else:
                            continue
                gscnv.close()
            #different row in filtered.xls and cnv.xls
            fileopen.close()
            filesave1.close()
            filesave2.close()
            os.system('sort -u '+self.out+'.'+chrom+'.sv.merged50.tmp.xls > '+self.out+'.'+chrom+'.sv.merged50.xls')
            os.system('sort -u '+self.out+'.'+chrom+'.sv.merged1050.tmp.xls > '+self.out+'.'+chrom+'.sv.merged1050.xls')
        print 'compare with gscnv set is done'
        
    def main(self):
        self.filterType()
        self.filterCount()
        self.filterRefgap()
        self.filterSegdup()
        self.filterHc()

        print 'lumpy filter is done at %s' % datetime.now()
if __name__ == '__main__':
    print 'start at %s' % datetime.now()
    Run = lumpyFilter()
    run = Run.main()

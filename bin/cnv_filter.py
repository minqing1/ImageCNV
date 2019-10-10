#! /usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
from datetime import datetime
def hundredFilter(samplename):
    c = 5
    path = 'example/'+samplename+'.merged.cnv.xls'
    sv_path = 'database/wgs_cnvfreq/sv_freq_V3.sorted.txt'
    filtered_path = 'example/'+samplename+'.merged.final.test.xls'
    fileopen = open(path,'r')
    filesave = open(filtered_path,'w')
    for i in fileopen:
        if i.startswith('chr'):
            linelist = i.split()
            sam_type = linelist[14].split(':')[1].lower()
            info = linelist[17].split(':')
            chrom = info[1]
            pos_start = int(info[2].split(';')[0])
            pos_end = int(info[3])
            sv = open(sv_path,'r')
            num = 0
            for j in sv:
                svlist = j.split()
                sv_chrom = svlist[0]
                sv_start = int(svlist[1].split('-')[0])
                sv_end = int(svlist[1].split('-')[1])
                sv_type = svlist[3]
                count = int(svlist[4])
                if chrom == sv_chrom and sam_type == sv_type:
                    if pos_start <= sv_start and sv_start<= pos_end <= sv_end:
                        overlap_length = pos_end - sv_start
                        sv_length = sv_end - sv_start
                        cnv_length = pos_end - pos_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5:
                            num += count
                    elif pos_end >= sv_end and sv_start <= pos_start <= sv_end:
                        overlap_length = sv_end - pos_start
                        sv_length = sv_end - sv_start
                        cnv_length = pos_end - pos_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5:
                            num += count
                    elif pos_start <= sv_start < sv_end <= pos_end:
                        overlap_length = sv_end - sv_start
                        cnv_length = pos_end - pos_start
                        sv_length = sv_end - sv_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5 :
                            num += count
                    elif pos_start == sv_start and pos_end == sv_end :
                        num += count
                    elif sv_start <= pos_start < pos_end <= sv_end:
                        overlap_length = pos_end - pos_start
                        sv_length = sv_end - sv_start
                        cnv_length = pos_end - pos_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5 :
                            num += count
                    else:
                        continue
                else:
                    continue
            if num <= c:
                filesave.write(i)
            sv.close()
        else:
            linelist = i.split()
            sam_type = linelist[2]
            info = linelist[3].split(':')
            chrom = info[0]
            pos_start = int(info[1].split('-')[0])
            pos_end = int(info[1].split('-')[1])
            sv = open(sv_path,'r')
            num = 0
            for j in sv:
                svlist = j.split()
                sv_chrom = svlist[0]
                sv_start = int(svlist[1].split('-')[0])
                sv_end = int(svlist[1].split('-')[1])
                sv_type = svlist[3]
                count = int(svlist[4])
                if chrom == sv_chrom and sam_type == sv_type:
                    if pos_start <= sv_start and sv_start<= pos_end <= sv_end:
                        overlap_length = pos_end - sv_start
                        sv_length = sv_end - sv_start
                        cnv_length = pos_end - pos_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5:
                            num += count
                    elif pos_end >= sv_end and sv_start <= pos_start <= sv_end:
                        overlap_length = sv_end - pos_start
                        sv_length = sv_end - sv_start
                        cnv_length = pos_end - pos_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5 :
                            num += count
                    elif pos_start <= sv_start < sv_end <= pos_end:
                        overlap_length = sv_end - sv_start
                        cnv_length = pos_end - pos_start
                        sv_length = sv_end - sv_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5 :
                            num += count
                    elif pos_start == sv_start and pos_end == sv_end :
                        num += count
                    elif sv_start <= pos_start < pos_end <= sv_end:
                        overlap_length = pos_end - pos_start
                        sv_length = sv_end - sv_start
                        cnv_length = pos_end - pos_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5 :
                            num += count
                    else:
                        continue
                else:
                    continue
            if num <= c:
                filesave.write(i)
            sv.close()
    fileopen.close()
    filesave.close()

    print '%s final cnv is done at %s ' %(samplename,datetime.now())
def main():
    if len(sys.argv)<2:
        print 'python 100k_filter.py samplename'
        sys.exit()
    else:
        #for i in open(sys.argv[1],'r'):
            #samplename = i.strip()
        print datetime.now()
        samplename = sys.argv[1]
        hundredFilter(samplename)

if __name__ == '__main__':
    run = main()

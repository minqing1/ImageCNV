#! /usr/bin/env python
# -*- coding:utf-8 -*-
# Description:cnvnator cnv set merge lumpy cnv set
import os
import sys
from datetime import datetime
class cnvMerged(object):
    def __init__(self,samplename):
        self.samplename = samplename
        self.cnv = 'example/'+samplename+'.cnv.final.xls'
        self.sv = 'example/'+samplename+'.sv.final.result.xls'
        self.path = 'example/'
        self.out = 'example/'+samplename

    def merged(self):
        ''' cnvnator cnv merge lumpy cnv set'''

        fcnv = open(self.cnv,'r')
        cnv50merged = self.out+'.cnv50merged.xls'
        sv50merged_tmp = self.out+'.sv50merged.tmp.xls'
        sv50merged = self.out+'.sv50merged.xls'
        filesave1 = open(cnv50merged,'w')
        filesave2 = open(sv50merged_tmp,'w')
        for i in fcnv:
            cnvlist = i.split()
            pos = cnvlist[3].split(':')[1]
            #print pos
            cnv_chr = cnvlist[3].split(':')[0]
            cnv_start = int(pos.split('-')[0])
            cnv_end = int(pos.split('-')[1])
            cnv_type = cnvlist[2]
            #print cnv_start,cnv_end
            fsv = open(self.sv,'r')
            for j in fsv:
                svlist = j.split()
                svinfo = svlist[17].split(':')
                sv_chr = svinfo[1]
                sv_start = int(svinfo[2].split(';')[0])
                sv_end = int(svinfo[3])
                sv_type = svlist[14].split(':')[1].lower()
                if cnv_chr == sv_chr and cnv_type == sv_type:
                    #1.back overlap in cnv >50%
                    if cnv_start <= sv_start< cnv_end <= sv_end:
                        overlap_length = cnv_end - sv_start
                        sv_length = sv_end - sv_start
                        cnv_length = cnv_end - cnv_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5:
                            filesave1.write(i)
                            filesave2.write(j)
                    #2.forehead overlap in cnv >50%
                    elif sv_start <= cnv_start < sv_end <= cnv_end:
                        overlap_length = sv_end - cnv_start
                        sv_length = sv_end - sv_start
                        cnv_length = cnv_end - cnv_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5:
                            filesave1.write(i)
                            filesave2.write(j)
                    #3. sv include in cnv 100%
                    elif cnv_start < sv_start < sv_end < cnv_end:
                        overlap_length = sv_end - sv_start
                        cnv_length = cnv_end - cnv_start
                        sv_length = sv_end - sv_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5:
                            filesave1.write(i)
                            filesave2.write(j)
                    #4 cnv == sv
                    elif cnv_start == sv_start and cnv_end == sv_end:
                        filesave1.write(i)
                        filesave2.write(j)
                    #5 cnv in sv
                    elif sv_start < cnv_start < cnv_end < sv_end:
                        overlap_length = cnv_end - cnv_start
                        sv_length = sv_end - sv_start
                        cnv_length = cnv_end - cnv_start
                        overlap_cnv = overlap_length / float(cnv_length)
                        overlap_sv = overlap_length / float(sv_length)
                        if overlap_cnv >= 0.5 and overlap_sv >= 0.5:
                            filesave1.write(i)
                            filesave2.write(j)
                    else:
                        continue
                else:
                    continue
            fsv.close()
        filesave1.close()
        filesave2.close()
        fcnv.close()
        os.system('sort -u '+sv50merged_tmp+' > '+sv50merged)
        os.system('grep -xvf '+cnv50merged+' '+self.cnv+' > '+self.out+'.cnvlt50.xls')
        os.system('grep -xvf '+sv50merged+' '+self.sv+' > '+self.out+'.svlt50.xls')
        os.system('cat '+self.out+'.cnvlt50.xls '+self.out+'.svlt50.xls '+sv50merged +' > '+self.out+'.merged.cnv.xls ')
        print 'cnv merge done'
    
    def compare_gscnv(self):
        ''' compare with NA12878 gscnv set '''
        merged = self.out+'.merged.xls'
        merged50gs = self.out+'.merged50gs.tmp.xls'
        merged1050gs = self.out+'.merged1050gs.tmp.xls'
        f = open(merged,'r')
        fsave1 = open(merged50gs,'w')
        fsave2 = open(merged1050gs,'w')
        for i in f:
            linelist = i.split()
            if linelist[0] == 'NA':
                info = linelist[3]
                cnv_chr = info.split(':')[0]
                postion = info.split(':')[1]
                cnv_start = int(postion.split('-')[0])
                cnv_end = int(postion.split('-')[1])
                cnv_type = linelist[2]
                gscnv = open(self.deletion,'r')

                for j in gscnv:
                    gaplist = j.split()
                    gap_chr = gaplist[0]
                    gap_start = int(gaplist[1])
                    gap_end = int(gaplist[2])
                    gap_type = 'deletion'

                    if cnv_chr == gap_chr and cnv_type == gap_type:
                        #1.back overlap in cnv >50%
                        if cnv_start <= gap_start < cnv_end <= gap_end:
                            overlap_length = cnv_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                            elif 0.1<=overlap_cnv and 0.1<= overlap_gap :
                                fsave2.write(i)

                        elif gap_start <= cnv_start < gap_end <= cnv_end:
                            overlap_length = gap_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                            elif 0.1 <= overlap_cnv and 0.1 <= overlap_gap:
                                fsave2.write(i)
                            
                    #3. gap #include in cnv 100%
                        elif cnv_start < gap_start < gap_end < cnv_end:
                            overlap_length = gap_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                            elif 0.1<=overlap_cnv and 0.1<= overlap_gap :
                                fsave2.write(i)
                            
                    #4 cnv == gap
                        elif cnv_start == gap_start and cnv_end == gap_end:
                            fsave1.write(i)
                    #5 cnv in gap
                        elif gap_start < cnv_start < cnv_end <gap_end:
                            overlap_length = cnv_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                            elif 0.1<=overlap_cnv and 0.1<=overlap_gap :
                                fsave2.write(i)
                            
                        else:
                            continue
                    else:
                        continue
                gscnv.close()
            else:
                cnvinfo = linelist[17].split(':')
                cnv_start = int(cnvinfo[2].split(';')[0])
                cnv_end = int(cnvinfo[3])
                cnv_chr = cnvinfo[1]
                gscnv = open(self.deletion,'r')
                cnv_type = linelist[14].split(':')[1].lower()
                for j in gscnv:
                    gaplist = j.split()
                    gap_chr = gaplist[0]
                    gap_start = int(gaplist[1])
                    gap_end = int(gaplist[2])
                    gap_type = 'deletion'

                    if cnv_chr == gap_chr and cnv_type == gap_type:
                        #1.back overlap in cnv >50%
                        if cnv_start <= gap_start< cnv_end <= gap_end:
                            overlap_length = cnv_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                            elif 0.1<=overlap_cnv and 0.1<= overlap_gap :
                                fsave2.write(i)
                            
                        elif gap_start <= cnv_start < gap_end <= cnv_end:
                            overlap_length = gap_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                            elif 0.1 <= overlap_cnv  and 0.1 <= overlap_gap:
                                fsave2.write(i)
                            
                        elif cnv_start < gap_start < gap_end < cnv_end:
                            overlap_length = gap_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                            elif 0.1<=overlap_cnv and 0.1<=overlap_gap :
                                fsave2.write(i)
                            
                    #4 cnv == gap
                        elif cnv_start == gap_start and cnv_end == gap_end:
                            fsave1.write(i)
                    #5 cnv in gap
                        elif gap_start < cnv_start < cnv_end <gap_end:
                            overlap_length = cnv_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                            elif 0.1<=overlap_cnv and 0.1<=overlap_gap :
                                fsave2.write(i)
                            
                        else:
                            continue
                    else:
                        continue
                gscnv.close()
        f.close()
        fsave1.close()
        fsave2.close()
        # remove same cnv call both in 50% and 10_50% set,respectively
        os.system('sort -u '+merged50gs+' > '+self.out+'.merged50gs.xls')

        print 'merge with NA12878 gold standard cnv set is done '
    def compare_sscnv(self):
        ''' using merged with gscnv lower than 10% cnv set 
            compare with NA12878 silver standard cnv set 
            rank merged ss set more than 50% and lower than 50%
        '''
        gscnvlt10 = self.out+'.merged.gscnv_lt10.xls'
        merged50ss = self.out+'.merged50ss.tmp.xls'
        f =open(gscnvlt10,'r')
        fsave1 = open(merged50ss,'w')
        for i in f:
            linelist =i.split()
            if linelist[0] == 'NA':
                info = linelist[3]
                cnv_chr = info.split(':')[0]
                postion = info.split(':')[1]
                cnv_start = int(postion.split('-')[0])
                cnv_end = int(postion.split('-')[1])
                sscnv = open(self.Sscnv,'r')
                for j in sscnv:
                    gaplist = j.split()
                    gap_chr = 'chr'+gaplist[0]
                    gap_start = int(gaplist[1])
                    gap_end = int(gaplist[2])
                    if cnv_chr == gap_chr:
                        if cnv_start <= gap_start< cnv_end <= gap_end:
                            overlap_length = cnv_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                        elif gap_start <= cnv_start < gap_end <= cnv_end:
                            overlap_length = gap_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                        elif cnv_start < gap_start < gap_end < cnv_end:
                            overlap_length = gap_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                    #4 cnv == gap
                        elif cnv_start == gap_start and cnv_end == gap_end:
                            fsave1.write(i)
                    #5 cnv in gap
                        elif gap_start < cnv_start < cnv_end <gap_end:
                            overlap_length = cnv_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                        else:
                            continue
                    else:
                        continue
                sscnv.close()
            else:
                cnvinfo = linelist[17].split(':')
                cnv_start = int(cnvinfo[2].split(';')[0])
                cnv_end = int(cnvinfo[3])
                cnv_chr = cnvinfo[1]
                sscnv = open(self.Sscnv,'r')
                for j in sscnv:
                    gaplist = j.split()
                    gap_chr = 'chr'+gaplist[0]
                    gap_start = int(gaplist[1])
                    gap_end = int(gaplist[2])
                    if cnv_chr == gap_chr:
                        if cnv_start <= gap_start< cnv_end <= gap_end:
                            overlap_length = cnv_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                        elif gap_start <= cnv_start < gap_end <= cnv_end:
                            overlap_length = gap_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                        elif cnv_start < gap_start < gap_end < cnv_end:
                            overlap_length = gap_end - gap_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                    #4 cnv == gap
                        elif cnv_start == gap_start and cnv_end == gap_end:
                            fsave1.write(i)
                    #5 cnv in gap
                        elif gap_start < cnv_start < cnv_end <gap_end:
                            overlap_length = cnv_end - cnv_start
                            cnv_length = cnv_end - cnv_start
                            gap_length = gap_end - gap_start
                            overlap_cnv = overlap_length / float(cnv_length)
                            overlap_gap = overlap_length / float(gap_length)
                            if overlap_cnv >= 0.5 and overlap_gap >= 0.5:
                                fsave1.write(i)
                        else:
                            continue
                    else:
                        continue
                sscnv.close()
        f.close()
        fsave1.close()
        os.system('sort -u '+merged50ss+' > '+self.out+'.merged50ss.xls')
        os.system('grep -xvf '+merged50ss+' '+gscnvlt10+' > '+self.out+'.No_overlap.xls')
        print 'merged with NA12878 silver standard cnv set is done'
    
    def main(self):
        self.merged()
        #self.compare_gscnv()
        #self.compare_sscnv()
        print 'merged cnv is done at %s ' % datetime.now()
if __name__ == '__main__':
    import sys
    usage = 'python merge_cnvsv.py binsize'
    if len(sys.argv)<2:
        print usage
        sys.exit()
    print 'start is done at %s ' % datetime.now()
    binsize = str(sys.argv[1])
    Run = cnvMerged(binsize)
    run = Run.main()

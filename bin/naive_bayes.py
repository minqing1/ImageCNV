#!/usr/bin/env python

import sys
import os
import glob
import random
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score

def filter_run(cnvdir, nordir, deldir, mutdir, filfile, binsize):
	stad_res = []
	pred_res = []
	filout = open(filfile, 'w')
	for m in range(1, 23):
		chrom = "chr"+str(m)
		cnv_list = glob.glob(cnvdir+"/"+chrom+"/*bed")
		mut_list = glob.glob(mutdir+"/*"+chrom+".snp.gff")
		het = {}
		last = 0
		for mutline in open(mut_list[0]):
			if mutline.startswith('#'):
				continue
			spl = mutline.split('\t')
			dep = [int(x) for x in spl[5].split(",")]
			(min_dep, max_dep) = (min(dep), max(dep))
			if 3 * min_dep > max_dep:
				if (int(spl[1]) - last) <= int(binsize):
					if last in het:
						het.pop(last)
				else:
					het[int(spl[1])] = 1
				last = int(spl[1])			
		for line in open(cnv_list[0]):
			line = line.rstrip()
			spl = line.split("\t")
			if 'dup' in spl[3]:
				filout.write(line+"\t0\n")
				continue
			het_num = 0
			length = 0
			for a in range(int(spl[1]), (int(spl[2]) + 1)):
				length += 1
				if a in het:
					het_num += 1

			length_short = length // 100
			length_short /= float(10) 
			stad_res.append(int(spl[-1]))
			delete_prob = prob_extract(deldir, length_short, het_num, 'del')
			normal_prob = prob_extract(nordir, length_short, het_num, 'nor')
			result_tmp = 1
			if (normal_prob == -1) or (delete_prob == -1):
				if (1000 * het_num) > length:
					result_tmp = 0
			else:
				if normal_prob > 2 * delete_prob:
					result_tmp = 0
			filout.write(line+"\t"+str(result_tmp)+"\n")
			pred_res.append(result_tmp)

	filout.close()	

	pre = precision_score(stad_res, pred_res)
	rec = recall_score(stad_res, pred_res)
	outfile = 'example/stat_result.xls'
	output = open(outfile, 'w')
	output.write("Precision\t"+str(pre)+"\n")
	output.write("Recall\t"+str(rec)+"\n")
	output.close()
		
def prob_extract(datadir, length_short, het_num, tp):
	prob = -1
	filename = datadir + "/" + str(length_short) + ".frq.xls"
	if os.path.exists(filename):
		for line in open(filename):
			if line.startswith("#"):
				if int(line.split("\t")[0].replace("#", "")) < 1000:
					break
			else:
				line = line.rstrip()
				spl = line.split("\t")
				if not line.startswith("Frequence"):
					continue
				if len(spl) < (het_num + 2):
					if 'del' in tp:
						prob = 0
					elif 'nor' in tp:
						prob = 1
				else:
					prob = float(spl[(het_num+1)])
	return prob

if __name__ == '__main__':
	cnv_dir = sys.argv[1]
	nor_dir = sys.argv[2]
	del_dir = sys.argv[3]
	mut_dir = sys.argv[4]
	fil_file = sys.argv[5]
	bin_size = int(sys.argv[6])
	filter_run(cnv_dir, nor_dir, del_dir, mut_dir, fil_file, bin_size)
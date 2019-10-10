#!/usr/bin/env python

import sys
import os

def image_identify_run(infile, picdir, modeldir, resdir, sample):
	if not os.path.exists(resdir+'image_tmp'):
		os.makedirs(resdir+'image_tmp')
	filfile = resdir + "/" + sample + ".image.fil.xls"
	output = open(filfile, "w")
	for line in open(infile):
		line = line.rstrip()
		spl = line.split("\t")
		picfile = picdir+"/"+sample+"."+spl[0]+"."+spl[1]+"-"+spl[2]+"."+spl[3]+".jpg"
		resfile = resdir+"/image_tmp/"+sample+"."+spl[0]+"."+spl[1]+"-"+spl[2]+"."+spl[3]+".xls"
		os.system("cd "+modeldir+" && python3 label_image.py "+picfile+ " >../../"+resfile)
		del_prob = 0
		dup_prob = 0
		for line2 in open(resfile):
			line2 = line2.rstrip()
			spl2 = line2.split(" ")
			if line2.startswith("del"):
				del_prob = float(spl2[3].replace(")", ""))
			elif line2.startswith("dup"):
				dup_prob = float(spl2[3].replace(")", ""))
		image_type = ""
		if dup_prob < del_prob:
			image_type = "deletion"
		else:
			image_type = "duplication"
		if image_type == spl[3]:
			output.write(line+"\t1\n")
		else:
			output.write(line+"\t0\n")
	output.close

if __name__ == '__main__':
	in_file = sys.argv[1]
	pic_dir = sys.argv[2]
	model_dir = sys.argv[3]
	res_dir = sys.argv[4]
	sample_name = sys.argv[5]
	image_identify_run(in_file, pic_dir, model_dir, res_dir, sample_name)
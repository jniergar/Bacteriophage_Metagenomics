#!/usr/bin/python
import sys, os, re, subprocess

# 2018-11-16 Jessmyn Niergarth
# Purpose: Get taxids from blastn and blastx output files.

# Usage:
usage = "Usage: python "+__file__+" [dir containing dir with sample outputs]"
if len(sys.argv[1:]) != 1:
	sys.exit("Incorrect number of command line arguments:\n"+usage)
dir_path = sys.argv[1]
out_accn = re.sub(r"\.py$", ".outAccs", __file__)
out_prot = re.sub(r"\.py$", ".outProts", __file__)

# Go through blastn & blastx results files:
dirs = os.listdir(dir_path)
prot_to_accn = {}
OUTA = open(out_accn, 'w')
OUTP = open(out_prot, 'w')
for one_dir in sorted(dirs):
	if not os.path.isdir(os.path.join(dir_path,one_dir)):
		dirs.remove(one_dir)
	else:
		blastn_path = os.path.join(dir_path,one_dir,one_dir+"_blastn_nt_acc.txt") #RN-1_blastn_nt_acc.txt
		blastx_path = os.path.join(dir_path,one_dir,one_dir+"_blastx_10e-3.txt")
		if os.path.isfile(blastn_path) and os.path.isfile(blastx_path):
# Get accessions from blastn file:
			INN = open(blastn_path)
			for line in INN:
				line_a = line.split("\t")
				OUTA.write(one_dir+"_"+line_a[0]+"\t"+line_a[4].split("|")[3]+"\n")
			INN.close()
# Get protein identifiers from blastx file:
			INX = open(blastx_path)
			for line in INX:
				line_a = line.split("\t")
				OUTP.write(one_dir+"_"+line_a[0]+"\t"+line_a[4]+"\n")
				prot_to_accn[line_a[4]] = 0
			INX.close()
		else:
			print("Couldn't find one/both of: "+blastn_path+", "+blastx_path)
			continue
OUTA.close()
OUTP.close()
print("Number of proteins: "+str(len(prot_to_accn)))

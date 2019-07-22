#!/usr/bin/python
import sys, os, re, subprocess

# 2018-11-16 Jessmyn Niergarth
# To get taxids from blastn and blastx output files

# Usage:
usage = "Usage: python "+__file__+" [dir containing dir with sample outputs]"
if len(sys.argv[1:]) != 1:
	sys.exit("Incorrect number of command line arguments:\n"+usage)
dir_path = sys.argv[1]
#outfile = re.sub(r"\.py$", ".out", __file__)
out_accn = re.sub(r"\.py$", ".outAccs", __file__)
out_prot = re.sub(r"\.py$", ".outProts", __file__)
out_prID = re.sub(r"\.py$", ".outProtsIDed", __file__)

# Go through blastn & blastx results files:
dirs = os.listdir(dir_path)
#cont_to_accn = {}
#cont_to_prot = {}
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
# Get accessions[?] from blastn file:
			INN = open(blastn_path)
			for line in INN:
				line_a = line.split("\t")
#				cont_to_accn[one_dir+"_"+line_a[0]] = line_a[1].split("|")[3]
				OUTA.write(one_dir+"_"+line_a[0]+"\t"+line_a[4].split("|")[3]+"\n")
			INN.close()
# Get protein identifiers from blastx file:
			INX = open(blastx_path)
			for line in INX:
				line_a = line.split("\t")
#				cont_to_prot[one_dir+"_"+line_a[0]] = line_a[4]
				OUTP.write(one_dir+"_"+line_a[0]+"\t"+line_a[4]+"\n")
				prot_to_accn[line_a[4]] = 0
			INX.close()
		else:
			print("Couldn't find one/both of: "+blastn_path+", "+blastx_path)
			continue
OUTA.close()
OUTP.close()
print("Number of proteins: "+str(len(prot_to_accn)))

# Entrez to get (genome?) accessions for each protein identifier

#esearch -db protein -query "YP_006987711.1" | elink -target nuccore | efetch -format acc


sys.exit() #DON'T need to do the following, because of file /home/jniergar/Documents/ncbi-blast-2.4.0+/tax/prot.accession2taxid !


# Find an accession number for each genome:

## TEMP:
#prot_to_accn = {}
#prot_to_accn["NP_040583.1"] = 0
#prot_to_accn["NP_040591.1"] = 0
#prot_to_accn["NP_046653.1"] = 0
#prot_to_accn["NP_046669.1"] = 0
#prot_to_accn["NP_579881.1"] = 0

OUTPID = open(out_prID, 'w')
for prot in sorted(prot_to_accn):
	cmd = 'esearch -db protein -query "'+prot+'" -silent | elink -target nuccore -silent | efetch -format acc -silent'
	result = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
	if len(result) > 0:
		result = result.rstrip()
		result_a = result.split("\n")
		print(prot+"\t"+result_a[0]) #is often, but not always, NC_*
#		prot_to_accn[prot] = result_a[0]
		OUTPID.write(prot+"\t"+result_a[0]+"\n")
	else:
		print("Second attempt on: "+prot)
		cmd1 = 'esearch -db protein -query "'+prot+'" -silent | efetch -format taxonomy | grep "taxname"'
		result1 = subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
#		print(result1)
		searchterm2 = re.sub(r'^.*?"(.*?)".*$',r'\1',result1)
#		print(searchterm2)
		cmd2 = 'esearch -db genome -query "'+searchterm2+'" | elink -target nuccore | efetch -format acc'
		result2 = subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
		if len(result2) > 0:
			result2_a = result2.split("\n")
			print(prot+"\t"+result2_a[0])
			OUTPID.write(prot+"\t"+result2_a[0]+"\n")
		else:
			print("\tDidn't work: "+prot)
#esearch -db protein -query "NP_579881" | efetch -format taxonomy | grep taxname
#        taxname "Human immunodeficiency virus 1",
#^output
#
#esearch -db genome -query "Human immunodeficiency virus 1" | elink -target nuccore | efetch -format acc
#^may have many, just grab first accn???
OUTPID.close()

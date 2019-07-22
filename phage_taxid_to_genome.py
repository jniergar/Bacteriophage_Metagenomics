#!/usr/bin/python
import re, os, sys, subprocess, time

# 2019-01-09 Jessmyn Niergarth
# Purpose: Use outputs of get_RPKMs_actuallyRPKMs.py to find and print out phage genomes from NCBI.

# Requires:
taxidfile = "get_RPKMs_actuallyRPKMs_taxids.out"
linfile = "get_RPKMs_actuallyRPKMs.out"
sortfile = "taxid_categorization_redux.txt" #file of virus hosts

outfile1 = re.sub(r"\.py", ".fasta", __file__)

# Get taxids:
all_taxids = []
INTAX = open(taxidfile)
for line in INTAX:
	if re.match(r"#", line): #skip first line
		continue
	else:
		all_taxids.append(line.split("\t")[0]) #first column contained taxids
INTAX.close()

# Get categorization (sorting) file:
choose_taxon = {}
taxon_record = {}
INSORT = open(sortfile)
for line in INSORT:
	if re.match(r"#", line): #skip first comment line (header)
		continue
	line_l = line.strip().split("\t")
	if line_l[0] in choose_taxon and choose_taxon[line_l[0]] != line_l[1]:
		sys.exit("Duplicate taxon with nonmatching host identifiers: "+line_l[0])
	else:
		choose_taxon[line_l[0]] = line_l[1]
		if line_l[1] == 'b':
			taxon_record[line_l[0]] = 0
INSORT.close()

# Only keep taxids corresponding to bacteriophages, according to lineage file:
vir_taxids = []
vir_lin = {} #virus lineages
counter = -1
INLIN = open(linfile)
for line in INLIN:
	if re.match(r"#", line):
		continue
	else:
		counter += 1
		matches = 0
		if re.search(r"Viruses\*superkingdom", line):
			for elem in line.split("\t")[0].split("/"):
				if elem in choose_taxon:
					matches += 1
					if choose_taxon[elem] == 'b':
						vir_taxids.append(all_taxids[counter])
						vir_lin[all_taxids[counter]] = line.split("\t")[0]
						taxon_record[elem] += 1
			if matches == 0:
				sys.exit("No subcategory of virus identified for line: "+line)
			elif matches <1:
				print("Multiple matches for line: "+line)
INLIN.close()
print("Number of taxids total: "+str(len(all_taxids))+"\nNumber of phage taxids: "+str(len(vir_taxids)))

# Intermediate checkin:
for onetax in taxon_record:
	if taxon_record[onetax] == 0:
		print("No taxids found for taxon:\t"+onetax)
	else:
		print("As expected, found taxids for taxon:\t"+onetax)

# Use esearch (from Entrez Direct: E-utilities) to find genomes corresponding to taxids and print to fasta file, with taxid at start of ID lines:
OUT1 = open(outfile1, 'w')
failed_taxids = []
for onevir in vir_taxids:
	print("Doing taxid:\t"+onevir+"\tlineage: "+vir_lin[onevir])
	cmd = 'esearch -db nuccore -query "txid'+onevir+' AND refseq[filter]" -silent | efetch -format fasta'
	result = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
	time.sleep(.2)
	if len(result) > 0:
		result_a = result.split('\n')
		for elem in result_a:
			elem = re.sub(r"^>", ">"+onevir+"-", elem)
			OUT1.write(elem+"\n")
	else:
		print("\tDidn't get sequence for: "+onevir)
		failed_taxids.append(onevir)

print("Re-attempting taxids that didn't work:")
double_failed = []
for onevir in failed_taxids:
	 print("Doing taxid:\t"+onevir+"\tlineage: "+vir_lin[onevir])
	cmd = 'esearch -db nuccore -query "txid'+onevir+' AND refseq[filter]" -silent | efetch -format fasta'
	result = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
	time.sleep(.2)
	if len(result) > 0:
		result_a = result.split('\n')
		for elem in result_a:
			elem = re.sub(r"^>", ">"+onevir+"-", elem)
			OUT1.write(elem+"\n")
	else:
		print("\tDidn't get sequence for: "+onevir)
		double_failed.append(onevir)
OUT1.close()

print("Twice-failed taxids: ")
print("\n".join(double_failed))

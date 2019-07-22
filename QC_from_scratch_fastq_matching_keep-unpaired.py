#!/usr/bin/python
import sys, os, re

# "QC_from_scratch_fastq_matching_keep-unpaired.py"

# 2017-05-15 Jessmyn Niergarth
# 2017-10-17 Jessmyn Niergarth: translated from Perl to Python
# Purpose: Performs quality control (QC) on a FASTQ file of reads, specifically:
#	Removes reads below [#] bp long or above [#] bp long
#	Removes reads with over [#] N's in total
#	Removes reads with over [#] N's in a row (continuous)
#	Removes homopolymers longer than [#]
# 'Remove' means fully removing seqname and seq of reads, but keeps unpaired reads that pass QC in a separate file.
# Note: names get shortened, split at the first space

# Usage:
usage = "Usage: python "+__file__+" [input FASTQ Fwd] [input FASTQ Rev] [output dir]"
if len(sys.argv[1:]) != 3:
	sys.exit("Incorrect number command line arguments:\n"+usage)

# Quality control parameters:
max = 301 #removes reads longer than this
min = 20 #removes reads shorter than this
Ntotl = 4 #removes reads that have this many Ns total or more
Ncont = 2 #removes reads that have this many continuous Ns or more
hp = 10 #removes reads with homopolymers this long or longer (Ntotl is redundant if $Ncont >= $Ntotl)
hpM = hp-1; #(necessary because of regex later)

# Setup
file1=sys.argv[1]
file2=sys.argv[2]
outdir=sys.argv[3]
outfile1_t = file1.split("/")[len(file1.split("/"))-1] #removes path
outfile2_t = file2.split("/")[len(file2.split("/"))-1]
if not os.path.isdir(outdir):
	sys.exit("Not a directory: "+outdir)
if not re.search(r"/$",outdir):
	outdir = outdir+"/"
outfile1 = outdir+re.sub(r"\.fastq",r"_QC.fastq",outfile1_t,flags=re.IGNORECASE)
outfile2 = outdir+re.sub(r"\.fastq",r"_QC.fastq",outfile2_t,flags=re.IGNORECASE)
outfile_unpaired = outdir+re.sub(r"R1_(.*?)\.fastq",r"RS_\1_QC.fastq",outfile1_t,flags=re.IGNORECASE)

if os.path.isfile(outfile_unpaired):
	os.system("rm -f "+outfile_unpaired)
#	sys.exit("Please delete file else will append to it: "+outfile_unpaired)

# Go through both files and find any seq that need to be removed in either
def fastqQC(onefile): #usage: fastqQC([name of file]); output: [dictionary of seqs to rm]
	# Setup:
	counter,seq_totl,seq_rm_Ntotl,seq_rm_Ncont,seq_rm_hp,seq_rm_short,seq_rm_long = (0,)*7 #all set to 0
	temp_seqname = None
	Ntotl_regex = r"^(.*N.*){"+str(Ntotl)+r"}$"
	Ncont_regex = r"N{"+str(Ncont)+r"}"
	hpM_regex = r"(.)\1{"+str(hpM)+r"}"
	seqs_rm = {}

	IN = open(onefile)
	print ("Reading in file:",onefile)
	for line in IN:
		line = line.strip('\n')
		line = line.strip('\cM')
		counter+=1
		if counter == 1:
			if not re.search('^@',line):
				sys.exit("Expected sequence name at: "+line)
			seq_totl+=1
			line_array = line.split(" ")
			temp_seqname = line_array[0]
			if not temp_seqname in seqs_rm:
				seqs_rm[temp_seqname] = 0; #default is to keep sequences
			else:
				sys.exit("Error; Duplicate sequence name: "+temp_seqname)
		elif counter == 2:
			if len(line) < min:
				seq_rm_short+=1
				seqs_rm[temp_seqname]+=1
				continue
			elif len(line) > max:
				seq_rm_long+=1
				seqs_rm[temp_seqname]+=1
				continue
			elif re.search(Ntotl_regex, line, re.IGNORECASE):
				seq_rm_Ntotl+=1
				seqs_rm[temp_seqname]+=1
				continue
			elif re.search(Ncont_regex, line, re.IGNORECASE):
				seq_rm_Ncont+=1
				seqs_rm[temp_seqname]+=1
				continue
			elif re.search(hpM_regex, line, re.IGNORECASE):
				seq_rm_hp+=1
				seqs_rm[temp_seqname]+=1
				continue
			else:
				if seqs_rm[temp_seqname] != 0:
					sys.exit("Disagreement about keeping sequence:"+temp_seqname)
		elif counter == 4:
			counter = 0;
			temp_seqname = None
		elif counter > 4:
			sys.exit("Problem: counter went over 4 at line:"+line)
	IN.close()

	print ("Number of sequences removed from:",onefile,"because:")
	print ("\t<",str(min),"bases in length:",str(seq_rm_short))
	print ("\t>",str(max),"bases in length:",str(seq_rm_long))
	print ("\t>=",str(Ntotl),"N\'s overall:",str(seq_rm_Ntotl))
	print ("\t>=",str(Ncont),"N\'s in a row:",str(seq_rm_Ncont))
	print ("\t>=",str(hp),"of homopolymer:",str(seq_rm_hp))

	total = seq_rm_short+seq_rm_short+seq_rm_Ntotl+seq_rm_Ncont+seq_rm_hp
	print ("Out of:",str(seq_totl),"removed:",total)

	return seqs_rm

seqs_rm1 = fastqQC(file1)
print ("\n",)
seqs_rm2 = fastqQC(file2) #keys will be seqnames, values will be 0 = keep, 1 = remove
for one_name in seqs_rm1:
	if not one_name in seqs_rm2:
		print ("Exists in file 1 but not file 2:"+one_name)
if not set(seqs_rm1) == set(seqs_rm2):
	sys.exit("Problem - sequence names not identical between FASTQ files")

# Go through dictionaries, comparing file origins of bad seqs
overall_keep = {} #keys will be seqnames, 0 = keep, 1 = remove
keeping, uniq_1, uniq_2, shared = (0,)*4
for one_name in seqs_rm1: #already checked that %seqs_rm1 and %seqs_rm2 had same keys
	if one_name in overall_keep:
		sys.exit("Duplicate sequence name: "+one_name)
	if seqs_rm1[one_name] == 0 and seqs_rm2[one_name] == 0:
		overall_keep[one_name] = 0;	keeping+=1;
	elif seqs_rm1[one_name] == 1 and seqs_rm2[one_name] == 0:
		overall_keep[one_name] = "u2";	uniq_1+=1; #bad seq uniquely in file1 i.e. FWD - 'u' for 'unpaired'
	elif seqs_rm1[one_name] == 0 and seqs_rm2[one_name] == 1:
		overall_keep[one_name] = "u1";	uniq_2+=1; #bad seq uniquely in file2 i.e. REV; bad in rev, so want fwd
	elif seqs_rm1[one_name] == 1 and seqs_rm2[one_name] == 1:
		overall_keep[one_name] = 1;	shared+=1
	else:
		txt = "Issue; Weird val[s] for: "+one_name+" Fwd: "+str(seqs_rm1[one_name])+"\tRev: "+str(seqs_rm2[one_name])
		sys.exit(txt)
print ("\nNumber seq keeping:\t"+str(keeping));		print ("Seq bad only in file 1:\t"+str(uniq_1))
print ("Seq bad only in file 2:\t"+str(uniq_2));	print ("Seq bad in both files:\t"+str(shared))

# Go through both files again, printing those seq that pass all QC for both input fastQs
def writingout(infile,outfile,u_file):
	counter = 0;
	print_boolean = 0;

	OUTUP = open(outfile_unpaired,"a")
	OUT = open(outfile,"w")
	IN = open(infile,"r")
	for line in IN:
		line = line.strip('\n')
		line = line.strip('\cM')
		counter+=1
		if counter == 1:
			if not re.search(r"^@", line):
				sys.exit("Expected sequence name at:"+line)
			line_array = line.split(" ")
			if overall_keep[line_array[0]] == 0:
				print_boolean = 1
				OUT.write(''.join([line,"\n"]))
			elif overall_keep[line_array[0]] == u_file:
				print_boolean = "u"
				OUTUP.write(''.join([line,"\n"]))
			else:
				print_boolean=0
		elif counter > 1 and counter < 5: #should match seq, '+', and quality lines
			if print_boolean == 1:
				OUT.write(''.join([line,"\n"]))
			elif print_boolean == "u":
				OUTUP.write(''.join([line,"\n"]))
			if counter == 4:
				counter = 0
	IN.close()
	OUT.close()
	OUTUP.close()

writingout(file1, outfile1, "u1") #(seqs bad only in file2 were marked as u1: meaning they are good in file1)
writingout(file2, outfile2, "u2")

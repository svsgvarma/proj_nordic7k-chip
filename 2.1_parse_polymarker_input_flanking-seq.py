#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""
Parse polymarker input and output a fasta file for blast.

./2.1_parse_polymarker_input_flanking-seq.py example_polymarker-1k.tsv /home/varma/proj/proj_nordic7k-chip/data/test_data/

"""

import sys

iupac = {"[A/G]": "R", "[G/A]": "R", "[C/T]": "Y", "[T/C]": "Y", "[G/C]": "S", "[C/G]": "S", "[A/T]": "W", "[T/A]": "W", "[G/T]": "K", "[T/G]": "K", "[A/C]": "M", "[C/A]": "M"}


def main():
	polymarker_input = sys.argv[1]
	workdir=sys.argv[2]
	outfile = polymarker_input.split(".")[0]
	#out = open(workdir+outfile+".fa", "w")
	outdown = open(workdir+outfile+"-down.fa", "w")
	outup = open(workdir+outfile+"-up.fa", "w")
	
	#### empty strings
	seqdownL50=0; seqdownL5=0; seqdownTot=0
	sequpL50=0; sequpL5=0; sequpTot=0
	
	# get snp position
	for line in open(workdir+polymarker_input):
		line = line.strip()
		if not line:
			continue
		snpname, seq = line.replace(" ","").split("\t")
		snpname = snpname.replace("_", "-") # in case there is already "_" in the snp name
		seq = seq.strip() # in case there is space in the input file
		pos = seq.find("[")
		snp = iupac[seq[pos:pos+5]]
		seq2 = seq[:pos] + snp + seq[pos+5:]
		seqdown = seq[:pos]
		seqdown50 = seqdown[-50:]
		sequp = seq[pos+5:]
		sequp50 = sequp[:50]

		if ((len(seqdown) < 50 ) and (len(seqdown) > 5)):
			#print(len(seqdown),len(snp),len(sequp) , str("Seqdown50:"), len(seqdown50))
			outdown.write(snpname + "_" + snp + "\n" + seqdown50 + "\n")
			seqdownL50+=1
		elif (len(seqdown) < 5):
		 	#print(len(seqdown),"too short DOWN read to BLAST")
		 	seqdownL5+=1
		else:
			#print(len(seqdown),len(snp),len(sequp) , str("Seqdown50:"), len(seqdown50))
			outdown.write(snpname + "_" + snp + "\n" + seqdown50 + "\n")
			seqdownTot+=1

		####
		if ((len(sequp) < 50 ) and (len(sequp) > 5)):
			#print(len(seqdown),len(snp),len(sequp) , str("Sequp50:"), len(sequp50))
			outup.write(snpname + "_" + snp + "\n" + sequp50 + "\n")
			sequpL50+=1
		elif (len(sequp) < 5):
		 	#print(len(sequp),"too short UP read to BLAST")
		 	sequpL5+=1
		else:
			#print(len(seqdown),len(snp),len(sequp) , str("Sequp50:"), len(sequp50))
			outup.write(snpname + "_" + snp + "\n" + sequp50 + "\n")
			sequpTot+=1

	#### stats about length 50 sequences
	print("seqdownL50: "+str(seqdownL50)+"\n"+ "seqdownL5: "+str( seqdownL5)+"\n"+"seqdownTot: "+str(seqdownTot)+"\n"+"sequpL50: "+str(sequpL50)+"\n"+"sequpL5: "+str(sequpL5)+"\n"+"sequpTot: "+str(sequpTot))
	outdown.close()
	outup.close()
	
	return 0

if __name__ == '__main__':
	main()

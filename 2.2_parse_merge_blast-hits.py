#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""
Parse polymarker input and output a fasta file for blast.

./2.2_parse_merge_blast-hits.py Nordic_Oat_SNP_sequences_AHA_230803_head1k.tsv /home/varma/proj/proj_nordic7k-chip/work/1.1_pars-BLAST/

"""

import sys
import re
import subprocess
from subprocess import *
from subprocess import call

iupac = {"[A/G]": "R", "[G/A]": "R", "[C/T]": "Y", "[T/C]": "Y", "[G/C]": "S", "[C/G]": "S", "[A/T]": "W", "[T/A]": "W", "[G/T]": "K", "[T/G]": "K", "[A/C]": "M", "[C/A]": "M"}

def main():
	polymarker_input = sys.argv[1]
	workdir=sys.argv[2]
	outfile = polymarker_input.split(".")[0]
	out = open(workdir+outfile+"_blast-merge.tsv", "w")

	#####--Function to search blast results---
	def srchdb_blastout(GName,FPATH):
		empty_str = str(".;.;.;.;.;.;.;.;.;.;.")
		try:
			True
			cmdFls1 = "grep -m 1 '"+str(GName)+"' "+str(FPATH)+""
			cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
			AnnoPars = cmdFls2.strip().decode().split("\t")
			# Replace "chr" with "str"
			for i in range(len(AnnoPars)):
				if (AnnoPars[i] == "chr"):
					AnnoPars[i] = empty_str
			grepout = str("\t".join(AnnoPars))
		except:
			False
			grepout = str("."+"\t"+empty_str+"\t"+empty_str+"\t"+empty_str)
		return grepout

	#####
	# write headers
	expfile_line1 = open(workdir+polymarker_input,'r').readline().strip()
	expfile_linesall = open(workdir+polymarker_input,'r').readlines()[1:]
	blastheader1=str(":(qseqid)")
	blastheader2=str(":(sseqid;pident;length;mismatch;gapopen;qstart;qend;sstart;send;evalue;bitscore")
	blastheader3=str(":(sseqid;pident)")
	outall_header = str(expfile_line1+"	Sequence-Down(-50)	Sequence-Up(+50)	Full-seq-len(seqdown;snppos;sequp)	Crop50-seq-len(seqdown50;snppos;sequp50)	All-seq_BLAST-hit1"+blastheader3+"	All-seq_BLAST-hit2"+blastheader3+"	All-seq_BLAST-hit3"+blastheader3+"	Seq-down_BLAST-hit1"+blastheader3+"	Seq-down_BLAST-hit2"+blastheader3+"	Seq-down_BLAST-hit3"+blastheader3+"	Seq-up_BLAST-hit1"+blastheader3+"	Seq-up_BLAST-hit2"+blastheader3+"	Seq-up_BLAST-hit3"+blastheader3+"	CODE_IUPAC_all-seq_BLAST-hit1"+blastheader1+"	all-seq_BLAST-hit1"+blastheader2+"	all-seq_BLAST-hit2"+blastheader2+"	all-seq_BLAST-hit3"+blastheader2+"	CODE_IUPAC_seq-down_BLAST-hit1"+blastheader1+"	seq-down_BLAST-hit1"+blastheader2+"	seq-down_BLAST-hit2"+blastheader2+"	seq-down_BLAST-hit3"+blastheader2+"	CODE_IUPAC_seq-up_BLAST-hit1"+blastheader1+"	seq-up_BLAST-hit1"+blastheader2+"	seq-up_BLAST-hit2"+blastheader2+"	seq-up_BLAST-hit3"+blastheader2+""+"\n")
	out.write(outall_header)
	# get blast hits position
	for line in expfile_linesall:
		lines = line.strip()
		lines_rmsps = "\t".join(lines.split())
		lines_snpid = lines.split()[0]
		##################
		# Count the number of DOWN and UP sequence lengths 
		if not line:
			continue
		seq = line.split("\t")[1]
		pos = seq.find("[")
		snp = iupac[seq[pos:pos+5]]
		seq2 = seq[:pos] + snp + seq[pos+5:]
		seqdown = seq[:pos]
		seqdown50 = seqdown[-50:]
		sequp = seq[pos+5:]
		sequp50 = sequp[:50]
		tol_seqlen = str(str(len(seqdown))+";"+str(len(snp))+";"+str(len(sequp)))
		crp50_seqlen = str(str(len(seqdown50))+";"+str(len(snp))+";"+str(len(sequp50)))

		##################
		# BLAST hit seqeunce order
		wdirblast= "/".join(workdir.split("/")[:-2])+"/1.0_annotate-BLAST/"
		order_seqmatch_all = srchdb_blastout(str(lines_snpid), wdirblast+"Nordic_Oat_SNP_sequences_AHA_230803_clusters.tsv")
		order_seqmatch_seqdown = srchdb_blastout(str(lines_snpid), wdirblast+"Nordic_Oat_SNP_sequences_AHA_230803_seq-down_clusters.tsv")
		order_seqmatch_sequp = srchdb_blastout(str(lines_snpid), wdirblast+"Nordic_Oat_SNP_sequences_AHA_230803_seq-up_clusters.tsv")

		#################
		# take only first two fiels from blast outfile6
		def parsblast2flds(elements_all):
			element_2fld_join = []
			for element in elements_all.split("\t"):
				if (element[0:5] != "AsSNP"):
					element_2fld = ";".join(element.split(";")[0:2])
					element_2fld_join.append(element_2fld)
			return str("\t".join(element_2fld_join))
		#################
		seqmatch_all_p2flds = parsblast2flds(order_seqmatch_all)
		seqmatch_seqdown_p2flds = parsblast2flds(order_seqmatch_seqdown) 
		seqmatch_sequp_p2flds = parsblast2flds(order_seqmatch_sequp) 

		#################
		outall_results= str(lines_rmsps+"\t"+seqdown50+"\t"+sequp50+"\t"+tol_seqlen+"\t"+crp50_seqlen+"\t"+seqmatch_all_p2flds+"\t"+seqmatch_seqdown_p2flds+"\t"+seqmatch_sequp_p2flds+"\t"+order_seqmatch_all+"\t"+ order_seqmatch_seqdown +"\t"+ order_seqmatch_sequp +"\n")
		##################
		# write output file
		out.write(outall_results)

	out.close()
	return 0

if __name__ == '__main__':
	main()

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
import numpy as np
import random
import pandas as pd
import itertools

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

	#####--Function to search 7K chip SNP data ---
	def srchdb_7k(GName,FPATH):
		try:
			True
			cmdFls1 = "grep -m 1 '"+str(GName)+"' "+str(FPATH)+""
			cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
			AnnoPars = cmdFls2.strip().decode().split("\t")
			# 
			grepout = str("\t".join(AnnoPars))
		except:
			False
			grepout = str(".	.	.	.	.")
		return grepout

	#####
	# write headers
	expfile_line1 = open(workdir+polymarker_input,'r').readline().strip()
	expfile_linesall = open(workdir+polymarker_input,'r').readlines()[1:]
	blastheader1=str(":(qseqid)")
	blastheader2=str(":(sseqid;pident;length;mismatch;gapopen;qstart;qend;sstart;send;evalue;bitscore")
	blastheader3=str(":(sseqid;pident;length)")
	exact_pos_join_header=str("Exact_pos_single(>=100pid)	Exact_pos_multi(>=100pid)	Exact_pos_single(<100pid)	Exact_pos_multi(<100pid)")
	sevenk_chip_header=str("SortOrder	CODE	Source	SNP_All	7K")
	outall_header = str(expfile_line1+"	Sequence-UP(-50)	Sequence-Down(+50)	Full-seq-len(sequp;snppos;seqdown)	Crop50-seq-len(sequp50;snppos;seqdown50)	All-seq_BLAST-hit1"+blastheader3+"	All-seq_BLAST-hit2"+blastheader3+"	All-seq_BLAST-hit3"+blastheader3+"	Seq-up_BLAST-hit1"+blastheader3+"	Seq-up_BLAST-hit2"+blastheader3+"	Seq-up_BLAST-hit3"+blastheader3+"	Seq-down_BLAST-hit1"+blastheader3+"	Seq-down_BLAST-hit2"+blastheader3+"	Seq-down_BLAST-hit3"+blastheader3+"	CODE_IUPAC_all-seq_BLAST-hit1"+blastheader1+"	all-seq_BLAST-hit1"+blastheader2+"	all-seq_BLAST-hit2"+blastheader2+"	all-seq_BLAST-hit3"+blastheader2+"	CODE_IUPAC_seq-up_BLAST-hits"+blastheader1+"	seq-up_BLAST-hit1"+blastheader2+"	seq-up_BLAST-hit2"+blastheader2+"	seq-up_BLAST-hit3"+blastheader2+"	CODE_IUPAC_seq-down_BLAST-hits"+blastheader1+"	seq-down_BLAST-hit1"+blastheader2+"	seq-down_BLAST-hit2"+blastheader2+"	seq-down_BLAST-hit3"+blastheader2+"")
	out.write(outall_header+"\t"+exact_pos_join_header+"\t"+sevenk_chip_header+"\n")

	totchrmatch=0; totchrunmatch=0
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
		
		sequp = seq[:pos]
		sequp50 = sequp[-50:]
		seqdown = seq[pos+5:]
		seqdown50 = seqdown[:50]

		tol_seqlen = str(str(len(sequp))+";"+str(len(snp))+";"+str(len(seqdown)))
		crp50_seqlen = str(str(len(sequp50))+";"+str(len(snp))+";"+str(len(seqdown50)))

		##################
		wdirblast= "/".join(workdir.split("/")[:-2])+"/1.0_annotate-BLAST/"
		order_seqmatch_all = srchdb_blastout(str(lines_snpid), wdirblast+"Nordic_Oat_SNP_sequences_AHA_230803_clusters.tsv")
		order_seqmatch_sequp = srchdb_blastout(str(lines_snpid), wdirblast+"Nordic_Oat_SNP_sequences_AHA_230803_seq-up_clusters.tsv")
		order_seqmatch_seqdown = srchdb_blastout(str(lines_snpid), wdirblast+"Nordic_Oat_SNP_sequences_AHA_230803_seq-down_clusters.tsv")
		#####
		order_srchdb_7k = srchdb_7k(str(lines_snpid), wdirblast+"Nordic_Oat_SNP_N7K_AHA_230907.tsv")

		#################
		# take only first two fields from blast outfile6
		def parsblast2flds(elements_all):
			element_2fld_join = []
			for element in elements_all.split("\t"):
				if (element[0:5] != "AsSNP") and (element[0:3] != "."):
					element_2fld = ";".join(element.split(";")[0:3])
					element_2fld_join.append(element_2fld)
			return str("\t".join(element_2fld_join))
		#################
		seqmatch_all_p2flds = parsblast2flds(order_seqmatch_all)
		seqmatch_sequp_p2flds = parsblast2flds(order_seqmatch_sequp)
		seqmatch_seqdown_p2flds = parsblast2flds(order_seqmatch_seqdown)

		#################
		# filter and summarized all fields from blast outfile6
		def parsblastflds(elements_all):
			element_flds_join = []
			for element in elements_all.split("\t"):
				if (element[0:5] != "AsSNP") and (element[0:3] != "."):
					element_flds = ";".join(element.split(";"))
					element_flds_join.append(element_flds)
			return str("\t".join(element_flds_join))
		#################

		#print(order_seqmatch_all)
		seqmatch_all_flds = parsblastflds(order_seqmatch_all)
		seqmatch_sequp_flds = parsblastflds(order_seqmatch_sequp)
		seqmatch_seqdown_flds = parsblastflds(order_seqmatch_seqdown)
		#print(seqmatch_all_flds)

		###############
		# filter and summarized 
		#print(seqmatch_all_p2flds,seqmatch_sequp_p2flds,seqmatch_seqdown_p2flds)
		#squp_pos_ = seqmatch_sequp_flds.split("\t")[0].split(";")[7:9]

		sqallchrlst =[]; sqallpidlst =[]; sqall_lenlst =[]; sqall_poslst=[]
		for elms in seqmatch_all_flds.split("\t"):
			sqallchrlst.append(elms.split(";")[0])
			sqallpidlst.append(elms.split(";")[1])
			sqall_lenlst.append(elms.split(";")[2])
			sqall_poslst.append(elms.split(";")[7:9])

		squpchrlst =[]; squppidlst =[]; squp_lenlst =[]; squp_poslst=[]
		for elms in seqmatch_sequp_flds.split("\t"):
			squpchrlst.append(elms.split(";")[0])
			squppidlst.append(elms.split(";")[1])
			squp_lenlst.append(elms.split(";")[2])
			squp_poslst.append(elms.split(";")[7:9])

		sqdownchrlst =[]; sqdownpidlst =[]; sqdown_lenlst =[]; sqdown_poslst=[]
		for elms in seqmatch_seqdown_flds.split("\t"):
			sqdownchrlst.append(elms.split(";")[0])
			sqdownpidlst.append(elms.split(";")[1])
			sqdown_lenlst.append(elms.split(";")[2])
			sqdown_poslst.append(elms.split(";")[7:9])

		#print(sqallchrlst[0],squpchrlst[0],sqdownchrlst[0])
		#print(sqallpidlst[0],squppidlst[0],sqdownpidlst[0])
		####################################################
		"""
		print("######################")
		print("###-first-UP-DOWN-####")
		print(squpchrlst[0],sqdownchrlst[0])
		print(squppidlst[0],sqdownpidlst[0])
		print(squp_lenlst[0],sqdown_lenlst[0])
		#
		print("###-second-UP-DOWN-####")
		print(squpchrlst[1],sqdownchrlst[1])
		print(squppidlst[1],sqdownpidlst[1])
		print(squp_lenlst[1],sqdown_lenlst[1])
		#
		print("###-third-UP-DOWN-####")
		print(squpchrlst[2],sqdownchrlst[2])
		print(squppidlst[2],sqdownpidlst[2])
		print(squp_lenlst[2],sqdown_lenlst[2])
		"""
		#####################
		# condition to replace no hits (.) to 0. 
		# Example: "100, 90, ."
		####
		# Function replace "." to 0.0 in pid / replace "." to 0 in length
		def replsdot(dist_list, str1, str2):
			dist_repdot = [float(i.replace(str1,str2)) if i == str1 else float(i) for i in dist_list]
			return dist_repdot

		# Function to sum all 100s in a list
		def sumall100(dist_pid):
			dist_sum = 0
			dist_pidsum = sum([dist_sum+1 if (i == 100) else dist_sum+ 0 for i in dist_pid])
			return dist_pidsum

		# Function to remove dot "." in a list
		def rmdot(listchr):
			ii = []
			[ii.append(i) if (i != ".") else i for i in listchr]
			return ii

		# Function to index "chr" string
		def index_string(list1, fstr):
			matched_indexes = []
			i = 0
			length = len(list1)
			while i < length:
				if fstr == list1[i]:
					matched_indexes.append(i)
				i += 1
			return matched_indexes

		# Function to Subject exact pos UP stream and range is reverse 
		#(example: ['116405498', '116405373'] and get len exact position)
		def pos_exact_UP(chr_idxlst,pos_idxlst, tol_seqlen):
			chr_pos_idxlst_sort =[]
			if (pos_idxlst[0] > pos_idxlst[1]):
				pos_idxlst_exact = str(int(pos_idxlst[1])-1)
				chr_pos_idxlst_sort.append(str(chr_idxlst)+":"+str(pos_idxlst_exact))
			else:
				pos_idxlst_exact = str(int(pos_idxlst[1])+1)
				chr_pos_idxlst_sort.append(str(chr_idxlst)+":"+str(pos_idxlst_exact))
			return chr_pos_idxlst_sort[0]

		# Function to Subject exact pos DOWN stream and range is reverse 
		#(example: ['116405498', '116405373'] and get len exact position)
		def pos_exact_DOWN(chr_idxlst,pos_idxlst, tol_seqlen):
			chr_pos_idxlst_sort =[]
			if (pos_idxlst[0] > pos_idxlst[1]):
				pos_idxlst_exact = str(int(pos_idxlst[0])+1)
				chr_pos_idxlst_sort.append(str(chr_idxlst)+":"+str(pos_idxlst_exact))
			else:
				pos_idxlst_exact = str(int(pos_idxlst[0])-1)
				chr_pos_idxlst_sort.append(str(chr_idxlst)+":"+str(pos_idxlst_exact))
			return chr_pos_idxlst_sort[0]

		#######################
		#print("###-UP-####")
		# Function to UP check and index all list
		def UP_allcond(up_pid_num, updown_chr, up_len_num):
			Exp_100pid_50len_lst_fun=[]; Exp_100pid_Less50len_lst_fun=[]
			for fstr in rmdot(updown_chr):
				up_chr_cln_idx = index_string(updown_chr, fstr)
				# check the pid is >100
				if (up_pid_num[up_chr_cln_idx[0]] == 100):
				# check the length is >50
					if (up_len_num[up_chr_cln_idx[0]] >= 50):
						squp_chr_idxlst = updown_chr[up_chr_cln_idx[0]]
						squp_pos_idxlst = squp_poslst[up_chr_cln_idx[0]]
						up_tol_seqlen = up_len_num[up_chr_cln_idx[0]]
						Exp_100pid_50len_lst_fun.append(pos_exact_UP(squp_chr_idxlst, squp_pos_idxlst, up_tol_seqlen))
					elif ( Exp_100pid_50len_lst_fun == []):
						if (up_len_num[up_chr_cln_idx[0]] < 50):
							squp_chr_idxlst = updown_chr[up_chr_cln_idx[0]]
							squp_pos_idxlst = squp_poslst[up_chr_cln_idx[0]]
							up_tol_seqlen = up_len_num[up_chr_cln_idx[0]]
							Exp_100pid_Less50len_lst_fun.append(pos_exact_UP(squp_chr_idxlst, squp_pos_idxlst, up_tol_seqlen))
			return Exp_100pid_50len_lst_fun,Exp_100pid_Less50len_lst_fun

		#######################
		#print("###-DOWN-####")
		# Function to DoWN check and index all list
		def DOWN_allcond(down_pid_num, updown_chr, down_len_num):
			Exp_100pid_50len_lst_fun=[]; Exp_100pid_Less50len_lst_fun=[]
			# Condition to check chr index	
			for fstr in rmdot(updown_chr):
				down_chr_cln_idx = index_string(updown_chr, fstr)
				# check the pid is >100
				#print(down_pid_num[down_chr_cln_idx[0]])
				if (down_pid_num[down_chr_cln_idx[0]] == 100):
					# check the length is >50
					if (down_len_num[down_chr_cln_idx[0]] >= 50):
						sqdown_chr_idxlst = updown_chr[down_chr_cln_idx[0]]
						sqdown_pos_idxlst = sqdown_poslst[down_chr_cln_idx[0]]
						down_tol_seqlen = down_len_num[down_chr_cln_idx[0]]
						Exp_100pid_50len_lst_fun.append(pos_exact_DOWN(sqdown_chr_idxlst,sqdown_pos_idxlst, down_tol_seqlen))
					elif ( Exp_100pid_50len_lst_fun == []):
						if (down_len_num[down_chr_cln_idx[0]] < 50):
							sqdown_chr_idxlst = updown_chr[down_chr_cln_idx[0]]
							sqdown_pos_idxlst = sqdown_poslst[down_chr_cln_idx[0]]
							down_tol_seqlen = down_len_num[down_chr_cln_idx[0]]
							Exp_100pid_Less50len_lst_fun.append(pos_exact_DOWN(sqdown_chr_idxlst,sqdown_pos_idxlst, down_tol_seqlen))
			return Exp_100pid_50len_lst_fun,Exp_100pid_Less50len_lst_fun
		#######################

		#######################
		#print("###-UP-####")
		# Function to UP check Less than 100 and index all list
		def UP_allcond_Less100(up_pid_num, updown_chr, up_len_num):
			Exp_Less100pid_50len_lst_fun=[]; Exp_Less100pid_Less50len_lst_fun=[]
			for fstr in rmdot(updown_chr):
				up_chr_cln_idx = index_string(updown_chr, fstr)
				# check the pid is >100
				if (up_pid_num[up_chr_cln_idx[0]] < 100 and up_pid_num[up_chr_cln_idx[0]] >= 95):
				# check the length is >50
					if (up_len_num[up_chr_cln_idx[0]] >= 50):
						squp_chr_idxlst = updown_chr[up_chr_cln_idx[0]]
						squp_pos_idxlst = squp_poslst[up_chr_cln_idx[0]]
						up_tol_seqlen = up_len_num[up_chr_cln_idx[0]]
						Exp_Less100pid_50len_lst_fun.append(pos_exact_UP(squp_chr_idxlst, squp_pos_idxlst, up_tol_seqlen))
					elif ( Exp_Less100pid_50len_lst_fun == []):
						if (up_len_num[up_chr_cln_idx[0]] < 50):
							squp_chr_idxlst = updown_chr[up_chr_cln_idx[0]]
							squp_pos_idxlst = squp_poslst[up_chr_cln_idx[0]]
							up_tol_seqlen = up_len_num[up_chr_cln_idx[0]]
							Exp_Less100pid_Less50len_lst_fun.append(pos_exact_UP(squp_chr_idxlst, squp_pos_idxlst, up_tol_seqlen))
			return Exp_Less100pid_50len_lst_fun,Exp_Less100pid_Less50len_lst_fun
		#######################
		#print("###-DOWN-####")
		# Function to DoWN check Less than 100 and index all list
		def DOWN_allcond_Less100(down_pid_num, updown_chr, down_len_num):
			Exp_Less100pid_50len_lst_fun=[]; Exp_Less100pid_Less50len_lst_fun=[]
			# Condition to check chr index	
			for fstr in rmdot(updown_chr):
				down_chr_cln_idx = index_string(updown_chr, fstr)
				# check the pid is >100
				#print(down_pid_num[down_chr_cln_idx[0]])
				if (down_pid_num[down_chr_cln_idx[0]] < 100 and down_pid_num[down_chr_cln_idx[0]] > 95):
					# check the length is >50
					if (down_len_num[down_chr_cln_idx[0]] >= 50):
						sqdown_chr_idxlst = updown_chr[down_chr_cln_idx[0]]
						sqdown_pos_idxlst = sqdown_poslst[down_chr_cln_idx[0]]
						down_tol_seqlen = down_len_num[down_chr_cln_idx[0]]
						Exp_Less100pid_50len_lst_fun.append(pos_exact_DOWN(sqdown_chr_idxlst,sqdown_pos_idxlst, down_tol_seqlen))
					elif ( Exp_Less100pid_50len_lst_fun == []):
						if (down_len_num[down_chr_cln_idx[0]] < 50):
							sqdown_chr_idxlst = updown_chr[down_chr_cln_idx[0]]
							sqdown_pos_idxlst = sqdown_poslst[down_chr_cln_idx[0]]
							down_tol_seqlen = down_len_num[down_chr_cln_idx[0]]
							Exp_Less100pid_Less50len_lst_fun.append(pos_exact_DOWN(sqdown_chr_idxlst,sqdown_pos_idxlst, down_tol_seqlen))
			return Exp_Less100pid_50len_lst_fun,Exp_Less100pid_Less50len_lst_fun

		# Pid number
		up_pid_num = replsdot(list([squppidlst[0],squppidlst[1],squppidlst[2]]),".","0.0")
		down_pid_num = replsdot(list([sqdownpidlst[0],sqdownpidlst[1],sqdownpidlst[2]]),".","0.0")
		# Pid number sum 
		updown_pidsum = list([sumall100(up_pid_num), sumall100(down_pid_num)])
		# Chr
		updown_chr = [list([squpchrlst[0],squpchrlst[1],squpchrlst[2]]),list([sqdownchrlst[0],sqdownchrlst[1],sqdownchrlst[2]])]
		updown_chr_join = [list([squpchrlst[0],squpchrlst[1],squpchrlst[2],sqdownchrlst[0],sqdownchrlst[1],sqdownchrlst[2]])]
		# Length number
		up_len_num = replsdot(list([squp_lenlst[0],squp_lenlst[1],squp_lenlst[2]]),".","0")
		down_len_num = replsdot(list([sqdown_lenlst[0],sqdown_lenlst[1],sqdown_lenlst[2]]),".","0")
		#######################################
		Exact_pos_100 =[]; Exact_pos_Less100 =[]
		########
		"""
		print("#######################")
		print(up_pid_num,down_pid_num)
		print("#### up  & down ")
		print(updown_chr[0],updown_chr[1])
		print(up_len_num, down_len_num)
		"""
		#Combinations with repetitions function
		if ((updown_pidsum[0] == 0) and (updown_pidsum[1] == 0)):
			True
			#print("no hits: up 0 & down 0!")
			if((sum(up_pid_num) != 0) and (sum(down_pid_num) == 0)):
				UP_allcond_lst = UP_allcond_Less100(up_pid_num, updown_chr[0], up_len_num)
				Exact_pos_Less100.append(list(set(UP_allcond_lst[0]) | set(UP_allcond_lst[1])))
			elif((sum(up_pid_num) == 0) and (sum(down_pid_num) != 0)):
				DOWN_allcond_lst = DOWN_allcond_Less100(down_pid_num, updown_chr[1], down_len_num)
				Exact_pos_Less100.append(list(set(DOWN_allcond_lst[0]) | set(DOWN_allcond_lst[1])))
			elif((sum(up_pid_num) == 0) and (sum(down_pid_num) == 0)):
				Exact_pos_Less100.append(["."])
			else:
				UP_allcond_Less100_lst = UP_allcond_Less100(up_pid_num, updown_chr[0], up_len_num)
				DOWN_allcond_Less100_lst = DOWN_allcond_Less100(down_pid_num, updown_chr[1], down_len_num)
				Exact_pos_Less100_50 = list(set(UP_allcond_Less100_lst[0]) | set(DOWN_allcond_Less100_lst[0]))
				Exact_pos_Less100_L50 = list(set(UP_allcond_Less100_lst[1]) | set(DOWN_allcond_Less100_lst[1]))
				if(Exact_pos_Less100_50 == []):
					Exact_pos_Less100.append(Exact_pos_Less100_L50)
				elif(Exact_pos_Less100_50 != []):
					Exact_pos_Less100.append(Exact_pos_Less100_50)
		elif ((updown_pidsum[0] == 0) and (updown_pidsum[1] == 1 or updown_pidsum[1] == 2 or updown_pidsum[1] == 3)):
			True
			# Condition to check chr index	
			DOWN_allcond_lst = DOWN_allcond(down_pid_num, updown_chr[1], down_len_num)
			Exact_pos_100.append(list(set(DOWN_allcond_lst[0]) | set(DOWN_allcond_lst[1])))

		elif ((updown_pidsum[0] == 1 or updown_pidsum[0] == 2 or updown_pidsum[0] == 3) and (updown_pidsum[1] == 0)):
			True
			# Condition to check chr index
			UP_allcond_lst = UP_allcond(up_pid_num, updown_chr[0], up_len_num)
			Exact_pos_100.append(list(set(UP_allcond_lst[0]) | set(UP_allcond_lst[1])))

		elif ((updown_pidsum[0] == 1 or updown_pidsum[0] == 2 or updown_pidsum[0] == 3) and (updown_pidsum[1] == 1 or updown_pidsum[1] == 2 or updown_pidsum[1] == 3)):
			True
			# Condition to check chr index
			UP_allcond_lst = UP_allcond(up_pid_num, updown_chr[0], up_len_num)
			DOWN_allcond_lst = DOWN_allcond(down_pid_num, updown_chr[1], down_len_num)

			UP_Exact_pos_100_50 = list(set(UP_allcond_lst[0]) | set(UP_allcond_lst[0]))
			UP_Exact_pos_100_L50 = list(set(UP_allcond_lst[1]) | set(UP_allcond_lst[1]))
			DOWN_Exact_pos_100_50 = list(set(DOWN_allcond_lst[0]) | set(DOWN_allcond_lst[0]))
			DOWN_Exact_pos_100_L50 = list(set(DOWN_allcond_lst[1]) | set(DOWN_allcond_lst[1]))
			#print(UP_Exact_pos_100_50, UP_Exact_pos_100_L50, DOWN_Exact_pos_100_50, DOWN_Exact_pos_100_L50)
			#print(UP_allcond_lst, DOWN_allcond_lst )
			if(UP_Exact_pos_100_50 == [] and DOWN_Exact_pos_100_50 == []):
				Exact_pos_100.append(list(set(UP_Exact_pos_100_L50) | set(DOWN_Exact_pos_100_L50)))
			elif(UP_Exact_pos_100_50 != [] and DOWN_Exact_pos_100_50 != []):
				Exact_pos_100.append(list(set(UP_Exact_pos_100_50) | set(DOWN_Exact_pos_100_50)))
			elif(UP_Exact_pos_100_50 != [] and DOWN_Exact_pos_100_50 == []):
				Exact_pos_100.append(UP_Exact_pos_100_50)
			elif(UP_Exact_pos_100_50 == [] and DOWN_Exact_pos_100_50 != []):
				Exact_pos_100.append(DOWN_Exact_pos_100_50)
		#################
		Exact_pos_100_single = []; Exact_pos_100_multi =[]
		Exact_pos_Less100_single = []; Exact_pos_Less100_multi =[]
		#print(Exact_pos_100,Exact_pos_Less100)
		#################
		#print(Exact_pos_100, up_pid_num, down_pid_num,  up_len_num, down_len_num)
		########-- 100 --########
		if(Exact_pos_100 != []):
			if(len(Exact_pos_100[0]) == 1):
				Exact_pos_100_single.append(Exact_pos_100[0][0])
			else:
				Exact_pos_100_single.append(".")
			###
			if(len(Exact_pos_100[0]) > 1):
				Exact_pos_100_multi.append(",".join(Exact_pos_100[0]))
			else:
				Exact_pos_100_multi.append(".")
		else:
			Exact_pos_100_single.append(".")
			Exact_pos_100_multi.append(".")
		#########-- Less100 --########
		if(Exact_pos_Less100 != []):
			if(len(Exact_pos_Less100[0]) == 1):
				Exact_pos_Less100_single.append(Exact_pos_Less100[0][0])
			else:
				Exact_pos_Less100_single.append(".")
			###
			if(len(Exact_pos_Less100[0]) > 1):
				Exact_pos_Less100_multi.append(",".join(Exact_pos_Less100[0]))
			else:
				Exact_pos_Less100_multi.append(".")
		else:
			Exact_pos_Less100_single.append(".")
			Exact_pos_Less100_multi.append(".")
		#################
		All_Exact_pos_join = str(Exact_pos_100_single[0])+"\t"+str(Exact_pos_100_multi[0])+"\t"+str(Exact_pos_Less100_single[0])+"\t"+str(Exact_pos_Less100_multi[0])
		#print(All_Exact_pos_join)
		##########################################################
		totchrunmatch +=1
		#################
		outall_results= str(lines_rmsps+"\t"+sequp50+"\t"+seqdown50+"\t"+tol_seqlen+"\t"+crp50_seqlen+"\t"+seqmatch_all_p2flds+"\t"+seqmatch_sequp_p2flds+"\t"+seqmatch_seqdown_p2flds+"\t"+order_seqmatch_all+"\t"+ order_seqmatch_sequp +"\t"+order_seqmatch_seqdown)
		#################
		# write output file
		out.write(outall_results+"\t"+All_Exact_pos_join+"\t"+order_srchdb_7k+"\n")
	out.close()
	#print("total-no matched to chr:"+str(totchrmatch)+"\n"+"total-no un-match to chr:"+str(totchrunmatch))
	print("script done!")
	return 0

if __name__ == '__main__':
	main()
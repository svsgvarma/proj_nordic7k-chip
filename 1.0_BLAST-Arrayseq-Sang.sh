#!/usr/bin/sh

#alpha version 0.1


############################
# Extract sequnces from the chip seq data

cat Nordic_Oat_SNP_sequences_AHA_230803.tsv | sed 1d | awk -F"\t" '{print ">"$1"\n"$2}' > Nordic_Oat_SNP_sequences_AHA_230803_seq.fa
cat Nordic_Oat_SNP_sequences_AHA_230803.tsv | sed 1d | awk -F"\t" '{print ">"$1"\t"$2}' > Nordic_Oat_SNP_sequences_AHA_230803_seq.tsv

#################
# Use pars ploymarker script 
#./parse_polymarker_input.py for_polymarker_1.tsv /home/varma/proj/proj_nordic7k-chip/data/test_data/
#./2.0_parse_polymarker_input.py Nordic_Oat_SNP_sequences_AHA_230803_seq.tsv /home/varma/proj/proj_nordic7k-chip/data/

############################
#echo "Annotaion with BLAST"

REF_PATH=/home/varma/proj/REF/Avena-sativa_Sangv1.1/
NCBIBLAST="/home/varma/softwares/BLAST/ncbi-blast-2.14.0+/bin"


#########################################################
REF=${REF_PATH}Asativa_sang_pseudomolecules
########
workdir="/home/varma/proj/proj_nordic7k-chip/"
indir=${workdir}data/
outdir=${workdir}work/1.0_annotate-BLAST/
mkdir -p $outdir
cd $outdir

#query=${indir}"Nordic_Oat_SNP_sequences_AHA_230803_seq.fa"
#outfl_1=${outdir}Nordic_Oat_SNP_sequences_AHA_230803_maxtar3_maxhsp1_eval1e05.tsv

query=${indir}"Nordic_Oat_SNP_sequences_AHA_230803_seq-up.fa"
outfl_1=${outdir}Nordic_Oat_SNP_sequences_AHA_230803_seq-up_maxtar3_maxhsp1_eval1e05.tsv

#query=${indir}"Nordic_Oat_SNP_sequences_AHA_230803_seq-down.fa"
#outfl_1=${outdir}Nordic_Oat_SNP_sequences_AHA_230803_seq-down_maxtar3_maxhsp1_eval1e05.tsv

########
#-perc_identity 95.0
#-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-05 

$NCBIBLAST/blastn -db $REF \
-query $query -num_threads 15 \
-max_target_seqs 3 -max_hsps 1 -evalue 1e-05 \
-outfmt 6 -out $outfl_1

echo "done BLAST script..."

#################
# Use pars ploymarker flaking sequence script

#./2.1_parse_polymarker_input_flanking-seq.py Nordic_Oat_SNP_sequences_AHA_230803_seq.tsv /home/varma/proj/proj_nordic7k-chip/data/


"""
# Some stats on Sequences
###--all sequences---
seqTot: 19600
seqTot blast hits: 19245

###-----
seqdownL50: 240
seqdownL5: 4837
seqdownG50: 14523
seqdownTot: 14763
seqdownTot blast hits: 14065

###-----
sequpL50: 227
sequpL5: 5945
sequpG50: 13428
sequpTot: 13655
sequpTot blast hits: 12969


"""
#####################

# Cluster based on first unique column and cluster all other columns keeping first unique column

#cat Nordic_Oat_SNP_sequences_AHA_230803_maxtar3_maxhsp3_eval1e05.tsv | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS=";" | awk '{split($12,a,";chr"); print $1"\tchr"a[2]"\tchr"a[3]"\tchr"a[4]}' FS="\t" OFS="\t" > Nordic_Oat_SNP_sequences_AHA_230803_clusters.tsv,
#cat Nordic_Oat_SNP_sequences_AHA_230803_seq-up_maxtar3_maxhsp1_eval1e05.tsv | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS=";" | awk '{split($12,a,";chr"); print $1"\tchr"a[2]"\tchr"a[3]"\tchr"a[4]}' FS="\t" OFS="\t" > Nordic_Oat_SNP_sequences_AHA_230803_seq-up_clusters.tsv
#cat Nordic_Oat_SNP_sequences_AHA_230803_seq-down_maxtar3_maxhsp1_eval1e05.tsv | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS=";" | awk '{split($12,a,";chr"); print $1"\tchr"a[2]"\tchr"a[3]"\tchr"a[4]}' FS="\t" OFS="\t" > Nordic_Oat_SNP_sequences_AHA_230803_seq-down_clusters.tsv

#####################

# Pars and merge all blast hits
# test
#time ./2.2_parse_merge_blast-hits.py Nordic_Oat_SNP_sequences_AHA_230803.tsv /home/varma/proj/proj_nordic7k-chip/work/1.1_pars-BLAST/
#
#real    4m57.054s
#user    1m50.555s
#sys     0m8.132s
#
###########
# Test run in the AzureVM
#time ./2.2_parse_merge_blast-hits.py Nordic_Oat_SNP_sequences_AHA_230803.tsv /home/bioinfoadmin/proj/proj_nordic7k-chip/work/1.1_pars-BLAST/
#
#real    4m23.512s
#user    2m5.809s
#sys     0m40.844s

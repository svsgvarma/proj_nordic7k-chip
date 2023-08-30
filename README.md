

## Workflow to extract sequences from the 7kchip data and BLAST against reference (Sang) genome data

### Pars Raw Arrayseq data 
./1.0_BLAST-Arrayseq-Sang.sh

### Parse ploymarker script is used to convert SNPs to IUPAC codes and add Up and Down streem flaking sequences
./2.1_parse_polymarker_input_flanking-seq.py Nordic_Oat_SNP_sequences_AHA_230803_seq.tsv /proj_nordic7k-chip/data/


### Pars and merge all blast hits
./2.2_parse_merge_blast-hits.py Nordic_Oat_SNP_sequences_AHA_230803.tsv /proj_nordic7k-chip/work/1.1_pars-BLAST/


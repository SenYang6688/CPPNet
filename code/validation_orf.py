import sys
import ORF
from Bio import SeqIO

def extract_feature_from_seq(seq, stt, stp):
	'''extract features of sequence from fasta entry'''
	
	stt_coden = stt.strip().split(',')
	stp_coden = stp.strip().split(',')
	# transtab = maketrans("ACGTNX","TGCANX")
	transtab = str.maketrans("ACGTNX", "TGCANX")
	mRNA_seq = seq.upper()
	mRNA_size = len(seq)
	tmp = ORF.ExtractORF(mRNA_seq)
	(CDS_size1, CDS_integrity, CDS_seq1) = tmp.longest_ORF(start=stt_coden, stop=stp_coden)
	return (mRNA_size, CDS_size1, CDS_integrity)


start_codons = 'ATG'
stop_codons = 'TAG,TAA,TGA'
Coverage = 0

seq_path = "/home/ys/work/cppred/code3_bak/CPPNet/data/validation_data/echi_all_ls_cp.fasta"
long_path = "./Echinococcus_granulosus_coding_RNA_validation.fa"
small_path = "./Echinococcus_granulosus_small_coding_RNA_validation.fa"

f_long = open(long_path,"w")
f_small = open(small_path,"w")

records = [record for record in SeqIO.parse(seq_path, "fasta")]

count_long = 0
count_small = 0
for i in range(len(records)):
	mRNA_size, CDS_size1, CDS_integrity = extract_feature_from_seq(records[i].seq, start_codons, stop_codons)
	if CDS_size1 > 303:
		count_long += 1
		f_long.write(">"+str(records[i].id)+"\n")
		f_long.write(str(records[i].seq)+"\n")
	else:
		count_small += 1
		f_small.write(">" + str(records[i].id) + "\n")
		f_small.write(str(records[i].seq) + "\n")
print(len(records))
print(count_long)
print(count_small)
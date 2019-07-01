import os, sys
import CTD
import ProtParam as PP
import ORF_length as my_len
from Bio import SeqIO
from Bio.Seq import Seq
import fickett
import FrameKmer
import argparse as agp
import rpy2.robjects as robjects
import pandas as pd
import numpy as np
from Bio.SeqUtils import ProtParam
from multiprocessing import Pool, cpu_count
import keras
from optparse import OptionParser
import warnings

warnings.filterwarnings("ignore")


class Fickett:
	'''
	calculate Fickett TESTCODE for full sequence
	NAR 10(17) 5303-531
	partially modified from source code of CPAT 1.2.1 downloaded from
	https://sourceforge.net/projects/rna-cpat/files/?source=navbar
	'''
	
	def __init__(self):
		self.content_parameter = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0]
		self.position_parameter = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
		'''
		newly calculated lookup table for RNA full length
		'''
		self.position_probability = {
			"A": [0.51, 0.55, 0.57, 0.52, 0.48, 0.58, 0.57, 0.54, 0.50, 0.36],
			"C": [0.29, 0.44, 0.55, 0.49, 0.52, 0.60, 0.60, 0.56, 0.51, 0.38],
			"G": [0.62, 0.67, 0.74, 0.65, 0.61, 0.62, 0.52, 0.41, 0.31, 0.17],
			"T": [0.51, 0.60, 0.69, 0.64, 0.62, 0.67, 0.58, 0.48, 0.39, 0.24],
		}
		self.position_weight = {"A": 0.062, "C": 0.093, "G": 0.205, "T": 0.154}
		self.content_probability = {
			"A": [0.40, 0.55, 0.58, 0.58, 0.52, 0.48, 0.45, 0.45, 0.38, 0.19],
			"C": [0.50, 0.63, 0.59, 0.50, 0.46, 0.45, 0.47, 0.56, 0.59, 0.33],
			"G": [0.21, 0.40, 0.47, 0.50, 0.52, 0.56, 0.57, 0.52, 0.44, 0.23],
			"T": [0.30, 0.49, 0.56, 0.53, 0.48, 0.48, 0.52, 0.57, 0.60, 0.51]
		}
		self.content_weight = {"A": 0.084, "C": 0.076, "G": 0.081, "T": 0.055}
	
	def look_up_position_probability(self, value, base):
		'''
		look up positional probability by base and value
		'''
		if float(value) < 0:
			return None
		for idx, val in enumerate(self.position_parameter):
			if (float(value) >= val):
				return float(self.position_probability[base][idx]) * float(self.position_weight[base])
	
	def look_up_content_probability(self, value, base):
		'''
		look up content probability by base and value
		'''
		if float(value) < 0:
			return None
		for idx, val in enumerate(self.content_parameter):
			if (float(value) >= val):
				return float(self.content_probability[base][idx]) * float(self.content_weight[base])
	
	def fickett_value(self, dna):
		'''
		calculate Fickett value from full RNA transcript sequence
		'''
		if len(dna) < 2:
			return 0
		fickett_score = 0
		dna = dna
		total_base = len(dna)
		A_content = float(dna.count("A")) / total_base
		C_content = float(dna.count("C")) / total_base
		G_content = float(dna.count("G")) / total_base
		T_content = float(dna.count("T")) / total_base
		
		phase_0 = dna[::3]
		phase_1 = dna[1::3]
		phase_2 = dna[2::3]
		
		phase_0_A = phase_0.count("A")
		phase_1_A = phase_1.count("A")
		phase_2_A = phase_2.count("A")
		phase_0_C = phase_0.count("C")
		phase_1_C = phase_1.count("C")
		phase_2_C = phase_2.count("C")
		phase_0_G = phase_0.count("G")
		phase_1_G = phase_1.count("G")
		phase_2_G = phase_2.count("G")
		phase_0_T = phase_0.count("T")
		phase_1_T = phase_1.count("T")
		phase_2_T = phase_2.count("T")
		
		A_content = float(phase_0_A + phase_1_A + phase_2_A) / total_base
		C_content = float(phase_0_C + phase_1_C + phase_2_C) / total_base
		G_content = float(phase_0_G + phase_1_G + phase_2_G) / total_base
		T_content = float(phase_0_T + phase_1_T + phase_2_T) / total_base
		A_position = np.max([phase_0_A, phase_1_A, phase_2_A]) / (np.min([phase_0_A, phase_1_A, phase_2_A]) + 1.0)
		C_position = np.max([phase_0_C, phase_1_C, phase_2_C]) / (np.min([phase_0_C, phase_1_C, phase_2_C]) + 1.0)
		G_position = np.max([phase_0_G, phase_1_G, phase_2_G]) / (np.min([phase_0_G, phase_1_G, phase_2_G]) + 1.0)
		T_position = np.max([phase_0_T, phase_1_T, phase_2_T]) / (np.min([phase_0_T, phase_1_T, phase_2_T]) + 1.0)
		
		fickett_score += self.look_up_content_probability(A_content, "A")
		fickett_score += self.look_up_content_probability(C_content, "C")
		fickett_score += self.look_up_content_probability(G_content, "G")
		fickett_score += self.look_up_content_probability(T_content, "T")
		
		fickett_score += self.look_up_position_probability(A_position, "A")
		fickett_score += self.look_up_position_probability(C_position, "C")
		fickett_score += self.look_up_position_probability(G_position, "G")
		fickett_score += self.look_up_position_probability(T_position, "T")
		
		return fickett_score


class FindCDS:
	'''
	Find the most like CDS in a given sequence
	The most like CDS is the longest ORF found in the sequence
	When having same length, the upstream ORF is printed
	modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar
	'''
	
	def __init__(self, seq):
		self.seq = seq
		self.result = (0, 0, 0, 0, 0)
		self.longest = 0
		self.basepair = {"A": "T", "T": "A", "U": "A", "C": "G", "G": "C", "N": "N", "X": "X"}
	
	def _reversecompliment(self):
		return "".join(self.basepair[base] for base in self.seq)[::-1]
	
	def get_codons(self, frame_number):
		'''
		Record every nucleotide triplet and its coordinate position for input sequence in one frame
		'''
		coordinate = frame_number
		while coordinate + 3 <= len(self.seq):
			yield (self.seq[coordinate:coordinate + 3], coordinate)
			coordinate += 3
	
	def find_longest_in_one(self, myframe, direction, start_codon, stop_codon):
		'''
		find the longest ORF in one reading myframe
		'''
		triplet_got = self.get_codons(myframe)
		starts = start_codon
		stops = stop_codon
		'''
		Extend sequence by triplet after start codon encountered
		End ORF extension when stop codon encountered
		'''
		while True:
			try:
				codon, index = triplet_got.__next__()
			except StopIteration:
				break
			if codon in starts and codon not in stops:
				'''
				find the ORF start
				'''
				orf_start = index
				end_extension = False
				while True:
					try:
						codon, index = triplet_got.__next__()
					except StopIteration:
						end_extension = True
						integrity = -1
					if codon in stops:
						integrity = 1
						end_extension = True
					if end_extension:
						orf_end = index + 3
						Length = (orf_end - orf_start)
						if Length > self.longest:
							self.longest = Length
							self.result = [direction, orf_start, orf_end, Length, integrity]
						if Length == self.longest and orf_start < self.result[1]:
							'''
							if ORFs have same length, return the one that if upstream
							'''
							self.result = [direction, orf_start, orf_end, Length, integrity]
						break
	
	def longest_orf(self, direction, start_codon={"ATG": None}, stop_codon={"TAG": None, "TAA": None, "TGA": None}):
		return_orf = ""
		for frame in range(3):
			self.find_longest_in_one(frame, "+", start_codon, stop_codon)
		return_orf = self.seq[self.result[1]:self.result[2]][:]
		start_coordinate = self.result[1]
		strand_direction = "+"
		orf_integrity = self.result[4]
		'''
		Also check reverse chain if -r is chosen
		'''
		if direction == "-":
			self.seq = self._reversecompliment()
			for frame in range(3):
				self.find_longest_in_one(frame, "-", start_codon, stop_codon)
			if self.result[0] == "-":
				return_orf = self.seq[self.result[1]:self.result[2]][:]
				start_coordinate = self.result[1]
				strand_direction = "-"
				orf_integrity = self.result[4]
		return return_orf, orf_integrity


def get_length(seq):
	return len(seq)


def get_fickett_value(seq):
	fickett = Fickett()
	return fickett.fickett_value(seq)


def get_stop_codon_num(seq):
	translate_prot = Seq(seq).translate()
	stop_num = translate_prot.count("*")
	return stop_num


def get_stop_codon_frequency(seq):
	stop_num = get_stop_codon_num(seq)
	transript_length = get_length(seq)
	stop_freq = float(stop_num) / transript_length
	return stop_freq


def get_orf(seq):
	findCDS = FindCDS(seq)
	return_orf, orf_integrity = findCDS.longest_orf(seq)
	# Whether the ORF starts with a START Codon and ends with a STOP Codon
	return return_orf, orf_integrity


def get_orf_coverge(seq):
	transript_length = get_length(seq)
	orf, _ = get_orf(seq)
	orf_length = len(orf)
	ORF_coverage = float(orf_length) / transript_length
	return ORF_coverage


def get_orf_frame_score(seq):
	ORF_length_in_frame1, _ = get_orf(seq)
	ORF_length_in_frame2, _ = get_orf(seq[1:])
	ORF_length_in_frame3, _ = get_orf(seq[2:])
	
	ORF_length_in_frame1 = len(ORF_length_in_frame1)
	ORF_length_in_frame2 = len(ORF_length_in_frame2)
	ORF_length_in_frame3 = len(ORF_length_in_frame3)
	
	ORF_len = [ORF_length_in_frame1, ORF_length_in_frame2, ORF_length_in_frame3]
	ORF_frame = ((ORF_len[0] - ORF_len[1]) ** 2 + (ORF_len[0] - ORF_len[2]) ** 2 + (ORF_len[1] - ORF_len[2]) ** 2) / 2
	return ORF_frame


def get_GC1(mRNA):
	if len(mRNA) < 3:
		numGC = 0
		mRNA = 'ATG'
	else:
		numGC = mRNA[0::3].count("C") + mRNA[0::3].count("G")
	return numGC * 1.0 / len(mRNA) * 3


def get_GC2(mRNA):
	if len(mRNA) < 3:
		numGC = 0
		mRNA = 'ATG'
	else:
		numGC = mRNA[1::3].count("C") + mRNA[1::3].count("G")
	return numGC * 1.0 / len(mRNA) * 3


def get_GC3(mRNA):
	if len(mRNA) < 3:
		numGC = 0
		mRNA = 'ATG'
	else:
		numGC = mRNA[2::3].count("C") + mRNA[2::3].count("G")
	return numGC * 1.0 / len(mRNA) * 3


def get_gc1_frame_score(seq):
	GC1_in_frame1 = get_GC1(seq)
	GC1_in_frame2 = get_GC1(seq[1:])
	GC1_in_frame3 = get_GC1(seq[2:])
	GC1_all = [GC1_in_frame1, GC1_in_frame2, GC1_in_frame3]
	GC1_frame = ((GC1_all[0] - GC1_all[1]) ** 2 + (GC1_all[0] - GC1_all[2]) ** 2 + (GC1_all[1] - GC1_all[2]) ** 2) / 2
	return GC1_frame


def get_gc2_frame_score(seq):
	GC2_in_frame1 = get_GC2(seq)
	GC2_in_frame2 = get_GC2(seq[1:])
	GC2_in_frame3 = get_GC2(seq[2:])
	GC2_all = [GC2_in_frame1, GC2_in_frame2, GC2_in_frame3]
	GC2_frame = ((GC2_all[0] - GC2_all[1]) ** 2 + (GC2_all[0] - GC2_all[2]) ** 2 + (GC2_all[1] - GC2_all[2]) ** 2) / 2
	return GC2_frame


def get_gc3_frame_score(seq):
	GC3_in_frame1 = get_GC3(seq)
	GC3_in_frame2 = get_GC3(seq[1:])
	GC3_in_frame3 = get_GC3(seq[2:])
	GC3_all = [GC3_in_frame1, GC3_in_frame2, GC3_in_frame3]
	GC3_frame = ((GC3_all[0] - GC3_all[1]) ** 2 + (GC3_all[0] - GC3_all[2]) ** 2 + (GC3_all[1] - GC3_all[2]) ** 2) / 2
	return GC3_frame


def get_stop_frame_score(seq):
	stop_num_in_frame1 = get_stop_codon_num(seq)
	stop_num_in_frame2 = get_stop_codon_num(seq[1:])
	stop_num_in_frame3 = get_stop_codon_num(seq[2:])
	stop_num_all = [stop_num_in_frame1, stop_num_in_frame2, stop_num_in_frame3]
	stop_num_frame = ((stop_num_all[0] - stop_num_all[1]) ** 2 + (stop_num_all[0] - stop_num_all[2]) ** 2 + (
			stop_num_all[1] - stop_num_all[2]) ** 2) / 2
	return stop_num_frame


def get_Mw(seq):
	translate_prot = ProtParam.ProteinAnalysis(
		str(Seq(seq).translate()).replace("*", "").replace("X", "").replace("Z", "").replace('B', '').replace('J', ''))
	mw = translate_prot.molecular_weight()
	return mw


def get_pI(seq):
	translate_prot = ProtParam.ProteinAnalysis(
		str(Seq(seq).translate()).replace("*", '').replace("X", "").replace("Z", "").replace('B', '').replace('J', ''))
	pI = translate_prot.isoelectric_point()
	return pI


def get_mw_div_pi(seq):
	mw = get_Mw(seq)
	pi = get_pI(seq)
	pi_mw = np.log10((float(mw) / pi) + 1)
	return pi_mw


def get_flexibility(seq):
	translate_prot = ProtParam.ProteinAnalysis(
		str(Seq(seq).translate()).replace("*", "").replace("X", "").replace("Z", "").replace('B', '').replace('J', ''))
	flexibility = translate_prot.flexibility()
	return flexibility


def get_gravy(seq):
	translate_prot = ProtParam.ProteinAnalysis(
		str(Seq(seq).translate()).replace("*", "").replace("X", "").replace("Z", "").replace('B', '').replace('J', ""))
	Gravy = translate_prot.gravy()
	return Gravy


def get_instablility_index(seq):
	translate_prot = ProtParam.ProteinAnalysis(
		str(Seq(seq).translate()).replace("*", "").replace("X", "").replace("Z", "").replace('B', '').replace('J', ''))
	instablility_index = translate_prot.instability_index()
	return instablility_index


def get_pi_mw_frame_score(seq):
	pi_mw_in_frame1 = get_mw_div_pi(seq)
	pi_mw_in_frame2 = get_mw_div_pi(seq[1:])
	pi_mw_in_frame3 = get_mw_div_pi(seq[2:])
	pi_mw_all = [pi_mw_in_frame1, pi_mw_in_frame2, pi_mw_in_frame3]
	pi_mw_frame = ((pi_mw_all[0] - pi_mw_all[1]) ** 2 + (pi_mw_all[0] - pi_mw_all[2]) ** 2 + (
			pi_mw_all[1] - pi_mw_all[2]) ** 2) / 2
	return pi_mw_frame


def feature_pipe(seq_path):
	res = []
	ids = []
	for seq in SeqIO.parse(seq_path, 'fasta'):
		ids.append(seq.id)
		seq = str(seq.seq)
		mw = get_Mw(seq)
		pi_mw = get_pi_mw_frame_score(seq)
		codon_num = get_stop_codon_num(seq)
		conon_f = get_stop_codon_frequency(seq)
		mv_div_pi = get_mw_div_pi(seq)
		gc1 = get_GC1(seq)
		gc2 = get_GC2(seq)
		gc3 = get_GC3(seq)
		gc1f = get_gc1_frame_score(seq)
		gc2f = get_gc2_frame_score(seq)
		gc3f = get_gc3_frame_score(seq)
		res.append(np.array([codon_num, conon_f, gc1, gc2, gc3, gc1f, gc2f, gc3f, pi_mw, mv_div_pi, mw]))
	return np.array(res), ids


def coding_nocoding_potential(input_file):
	coding = {}
	noncoding = {}
	for line in open(input_file).readlines():
		fields = line.split()
		if fields[0] == 'hexamer': continue
		coding[fields[0]] = float(fields[1])
		noncoding[fields[0]] = float(fields[2])
	return coding, noncoding


def output_feature(seq_file, hex_file):
	res = []
	coding, noncoding = coding_nocoding_potential(hex_file)
	for seq in SeqIO.parse(seq_file, 'fasta'):
		A, T, G, C, AT, AG, AC, TG, TC, GC, A0, A1, A2, A3, A4, T0, T1, T2, T3, T4, G0, G1, G2, G3, G4, C0, C1, C2, C3, C4 = CTD.CTD(
			seq.seq)
		insta_fe, PI_fe, gra_fe = PP.param(seq.seq)
		fickett_fe = fickett.fickett_value(seq.seq)
		hexamer = FrameKmer.kmer_ratio(seq.seq, 6, 3, coding, noncoding)
		Len, Cov, inte_fe = my_len.len_cov(seq.seq)
		tem = [A, T, G, C, AT, AG, AC, TG, TC, GC, A0, A1, A2, A3, A4, T0, T1, T2, T3, T4, G0, G1, G2, G3, G4, C0, C1,
		       C2, C3, C4, inte_fe, Cov, insta_fe, PI_fe, gra_fe, hexamer, fickett_fe, Len]
		res.append(tem)
	return np.array(res)


def get_r_feature(seq_path):
	r_code = '''
		library(LncFinder)
		library(seqinr)
		seq = read.fasta(file = {path})
		features = extract_features(seq, label = NULL, SS.features = FALSE,format = "DNA", frequencies.file = "human", parallel.cores = 8)
		features
		'''.format(path="'" + seq_path + "'")
	res = robjects.r(r_code)
	colnames = ["ORF.Max.Len",
	            "ORF.Max.Cov",
	            "Seq.lnc.Dist",
	            "Seq.pct.Dist",
	            "Seq.Dist.Ratio",
	            "Signal.Peak",
	            "SNR",
	            "Signal.Min",
	            "Signal.Q1",
	            "Signal.Q2",
	            "Signal.Max"]
	df = pd.DataFrame(np.array(res).T, columns=colnames).ix[:, 2:]
	return df


def get_feature(seq_path, hex_file):
	my_feature, ids = feature_pipe(seq_path)
	lncfinder_feature = get_r_feature(seq_path).values
	cpp_feature = output_feature(seq_path, hex_file)
	res = np.concatenate([cpp_feature, lncfinder_feature[:, 2:], my_feature], axis=1)
	features_names = ["F" + str(num) for num in range(1, 57, 1)]
	return res, ids, features_names


if __name__ == "__main__":
	hex_file = "./code/Cross_species_Hexamer.tsv"
	model_path = "./code/Cross_species_CPPNet.h5"
	
	parser = OptionParser()
	parser.add_option("-i", "--input_file", dest="input_file", help="Input fasta file path with fasta format")
	parser.add_option("-o", "--output_file", dest="output_file", help="Output file path with csv format")
	option, args = parser.parse_args()
	features, ids, features_names = get_feature(option.input_file, hex_file)
	features_names = ["protein-coding RNA Probability"] + features_names
	model = keras.models.load_model(model_path, compile=False)
	prob = model.predict(features)[:, 1].reshape((-1, 1))
	res = np.concatenate([prob, features], axis=1)
	df = pd.DataFrame(res, index=ids, columns=features_names)
	df.to_csv(option.output_file)
	print(df)




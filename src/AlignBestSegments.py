import os

import pandas as pd


CODON_LIST = [('A', 'GCT'), ('A', 'GCC'), ('A', 'GCA'), ('A', 'GCG'),
 ('L', 'TTA'), ('L', 'TTG'), ('L', 'CTT'), ('L', 'CTC'), ('L', 'CTA'), ('L', 'CTG'),
 ('R', 'CGT'), ('R', 'CGC'), ('R', 'CGA'), ('R', 'CGG'), ('R', 'AGA'), ('R', 'AGG'),
 ('K', 'AAA'), ('K', 'AAG'),
 ('N', 'AAT'), ('N', 'AAC'),
 ('M', 'ATG'),
 ('D', 'GAT'), ('D', 'GAC'),
 ('F', 'TTT'), ('F', 'TTC'),
 ('C', 'TGT'), ('C', 'TGC'),
 ('P', 'CCT'), ('P', 'CCC'), ('P', 'CCA'), ('P', 'CCG'),
 ('Q', 'CAA'), ('Q', 'CAG'),
 ('S', 'TCT'), ('S', 'TCC'), ('S', 'TCA'), ('S', 'TCG'), ('S', 'AGT'), ('S', 'AGC'),
 ('E', 'GAA'), ('E', 'GAG'),
 ('T', 'ACT'), ('T', 'ACC'), ('T', 'ACA'), ('T', 'ACG'),
 ('G', 'GGT'), ('G', 'GGC'), ('G', 'GGA'), ('G', 'GGG'),
 ('W', 'TGG'),
 ('H', 'CAT'), ('H', 'CAC'),
 ('Y', 'TAT'), ('Y', 'TAC'),
 ('I', 'ATT'), ('I', 'ATC'), ('I', 'ATA'),
 ('V', 'GTT'), ('V', 'GTC'), ('V', 'GTA'), ('V', 'GTG'),
 ('*', 'TAA'), ('*', 'TGA'), ('*', 'TAG')]


CODONS = {}
for aa, codon in CODON_LIST:
	CODONS[aa] = CODONS.get(aa, []) + [codon]


def align_nuc_to_aa(seq, gene):

	def _align_slice(codon_pos, gene_symbol, aminoacid, prev_aligned = [True] * 6):
		aligned = [False] * 6
		for i, codon in enumerate(CODONS[aminoacid]):
			aligned[i] = prev_aligned[i] & (codon[codon_pos] == gene_symbol)
		return sum(aligned) > 0, aligned

	score = 0

	flag = False

	for aa_pos, aa in enumerate(seq):
		if aa_pos*3 < len(gene):
			flag, aligned = _align_slice(0, gene[aa_pos*3], aa)

			if flag:
				score += 1

				if (aa_pos*3 + 1) < len(gene):
					flag, aligned = _align_slice(1, gene[aa_pos*3 + 1], aa, aligned)

					if flag:
						score += 1

						if (aa_pos*3 + 2) < len(gene):
							flag, aligned = _align_slice(2, gene[aa_pos*3 + 2], aa, aligned)

							score += flag
		if not flag:
			break

	return score


def align_nuc_to_aa_rev(seq, gene):
	
	def _align_slice(codon_pos, gene_symbol, aminoacid, prev_aligned = [True] * 6):
		aligned = [False] * 6
		for i, codon in enumerate(CODONS[aminoacid]):
			aligned[i] = prev_aligned[i] & (codon[codon_pos] == gene_symbol)
		return sum(aligned) > 0

	score = 0

	for aa_pos, aa in enumerate(reversed(seq)):
		if aa_pos*3 < len(gene):
			flag, aligned = _align_slice(2, gene[-(aa_pos*3 + 1)], aa)

			if flag:
				score += 1

				if (aa_pos*3 + 1) < len(gene):
					flag, aligned = _align_slice(1, gene[-(aa_pos*3 + 2)], aa, aligned)

					if flag:
						score += 1

						if (aa_pos*3 + 2) < len(gene):
							flag, aligned = _align_slice(0, gene[-(aa_pos*3 + 3)], aa, aligned)

							score += flag
		if not flag:
			break

	return score


def align_segments_and_write(sequence_files, segments_filepath="./segments.txt"):

	def _fix_v(index, row, df, gene_type, segments):
		seg_gene_type = ""
		if gene_type == "alpha":
			seg_gene_type = "TRA"
		elif gene_type == "beta":
			seg_gene_type = "TRB"
		else:
			print("Error: unknown gene type", gene_type)
			return 0

		fixed_seg = ""
		max_score = -1
		old_score = 0
		if not pd.isnull(row["cdr3." + gene_type]):
			for _, seg_row in segments[(segments.species == row["species"]) & (segments.gene == seg_gene_type) & (segments.segment == "Variable")].iterrows():
				cur_score = align_nuc_to_aa(row["cdr3." + gene_type], seg_row["seq"][seg_row["ref"] - 3:])
				if cur_score > max_score:
					max_score = cur_score
					fixed_seg = seg_row["id"]

				# TODO: search for a row after all this iterations. Same for J.
				if seg_row["id"][:seg_row["id"].find("*")] == row["v." + gene_type]:
					old_score = cur_score

			df.set_value(index, "v." + gene_type + ".fixed", fixed_seg)
			df.set_value(index, "v." + gene_type + ".end.fixed", max_score)
			df.set_value(index, "v." + gene_type + ".end.old", old_score)


	def _fix_j(index, row, df, gene_type, segments):
		seg_gene_type = ""
		if gene_type == "alpha":
			seg_gene_type = "TRA"
		elif gene_type == "beta":
			seg_gene_type = "TRB"
		else:
			print("Error: unknown gene type", gene_type)
			return 0

		fixed_seg = ""
		max_score = -1
		old_score = 0
		if not pd.isnull(row["cdr3." + gene_type]):
			for _, seg_row in segments[(segments.species == row["species"]) & (segments.gene == seg_gene_type) & (segments.segment == "Joining")].iterrows():
				cur_score = align_nuc_to_aa(row["cdr3." + gene_type], seg_row["seq"][:seg_row["ref"] + 4])
				if cur_score > max_score:
					max_score = cur_score
					fixed_seg = seg_row["id"]

				if seg_row["id"][:seg_row["id"].find("*")] == row["j." + gene_type]:
					old_score = cur_score

			df.set_value(index, "j." + gene_type + ".fixed", fixed_seg)
			df.set_value(index, "j." + gene_type + ".start.fixed", len(row["cdr3." + gene_type]) - max_score)
			df.set_value(index, "j." + gene_type + ".start.old", len(row["cdr3." + gene_type]) - old_score)


	if sequence_files is not list:
		if os.path.isdir(sequence_files):
			sequence_files = [sequence_files + "/" + x for x in os.listdir(sequence_files) if os.path.isfile(sequence_files + "/" + x)]
		else:
			sequence_files = [sequence_files]

	segments = pd.read_csv(segments_filepath, sep="\t")
	segments.columns = ["species", "gene", "segment", "id", "ref", "seq"]

	for file_i, seq_file in enumerate(sequence_files):
		print(file_i, "/", len(sequence_files), "Processing", seq_file, "...", end = "\t")
		df = pd.read_csv(seq_file, sep="\t")

		df["v.alpha.fixed"] = "NA"
		df["j.alpha.fixed"] = "NA"
		df["v.beta.fixed"] = "NA"
		df["j.beta.fixed"] = "NA"

		df["v.alpha.end.fixed"] = "-1"
		df["j.alpha.start.fixed"] = "-1"
		df["v.beta.end.fixed"] = "-1"
		df["j.beta.start.fixed"] = "-1"

		df["v.alpha.end.old"] = "-1"
		df["j.alpha.start.old"] = "-1"
		df["v.beta.end.old"] = "-1"
		df["j.beta.start.old"] = "-1"

		for index, row in df.iterrows():
			_fix_v(index, row, df, "alpha", segments)
			_fix_j(index, row, df, "alpha", segments)
			_fix_v(index, row, df, "beta", segments)
			_fix_j(index, row, df, "beta", segments)

		df.to_csv(seq_file[:-3] + "fixed.txt", sep="\t", index=False)
		print("Done.")


if __name__ == "__main__":
	print(align_nuc_to_aa("L", "TTC"))
	# align_segments_and_write("../chunks/")
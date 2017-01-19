import os
import sys
import json

import pandas as pd

from __future__ import print_function


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


def align_segments_and_write(full_table, segments_filepath="./segments.txt"):

	def _fix_json(index, row, df, gene_type, segments):
		seg_gene_type = ""
		if gene_type == "alpha":
			seg_gene_type = "TRA"
		elif gene_type == "beta":
			seg_gene_type = "TRB"
		else:
			print("Error: unknown gene type", gene_type)
			return 0

		if not pd.isnull(row["cdr3fix." + gene_type]):
			json_val = json.loads(row["cdr3fix." + gene_type])

			# если скоры одинаковые??

			# VARIABLE
			max_score = -1
			old_score = 0
			fixed_seg = "None"
			for _, seg_row in segments[(segments.species == row["species"]) & (segments.gene == seg_gene_type) & (segments.segment == "Variable")].iterrows():
				cur_score = align_nuc_to_aa(row["cdr3." + gene_type], seg_row["seq"][seg_row["ref"] - 3:])
				if cur_score > max_score:
					max_score = cur_score
					fixed_seg = seg_row["id"]

				# TODO: search for a row after all this iterations. Same for J.
				if seg_row["id"][:seg_row["id"].find("*")] == row["v." + gene_type]:
					old_score = cur_score

			# NoFixNeeded -> NoFix or ChangeSegment
			# FixAdd, FixReplace, FixTrim -> ChangedSequence
			# fail -> Failed

			# кроме тех где NoFixNeeded
			json_val["oldVId"] = json_val["vId"]
			json_val["vId"] = fixed_seg
			json_val["oldVEnd"] = json_val["vEnd"]
			json_val["vEnd"] = max_score
			json_val["oldjFixType"] = json_val["jFixType"]
			if max_score != -1:
				if json_val["jFixType"] == "NoFixNeeded":
					if json_val["vId"] == json_val["oldVId"]:
						json_val["jFixType"] = "NoFix"
					else:
						json_val["jFixType"] = "ChangeSegment"
				elif json_val["jFixType"] in ["FixAdd", "FixReplace", "FixTrim"]:
					json_val["jFixType"] = "ChangeSequence"
				else:
					json_val["jFixType"] = "Failed"
			else:
				json_val["jFixType"] = "Failed"

			# JOINING
			max_score = -1
			old_score = 0
			fixed_seg = "None"
			for _, seg_row in segments[(segments.species == row["species"]) & (segments.gene == seg_gene_type) & (segments.segment == "Joining")].iterrows():
				cur_score = align_nuc_to_aa(row["cdr3." + gene_type], seg_row["seq"][:seg_row["ref"] + 4])
				if cur_score > max_score:
					max_score = cur_score
					fixed_seg = seg_row["id"]

				if seg_row["id"][:seg_row["id"].find("*")] == row["j." + gene_type]:
					old_score = cur_score

			json_val["oldJId"] = json_val["jId"]
			json_val["jId"] = fixed_seg
			json_val["oldJStart"] = json_val["jStart"]
			json_val["jStart"] = len(row["cdr3." + gene_type]) - max_score
			json_val["oldJFixType"] = json_val["jFixType"]
			if max_score != -1:
				if json_val["jFixType"] == "NoFixNeeded":
					if json_val["vId"] == json_val["oldVId"]:
						json_val["jFixType"] = "NoFix"
					else:
						json_val["jFixType"] = "ChangeSegment"
				elif json_val["jFixType"] in ["FixAdd", "FixReplace", "FixTrim"]:
					json_val["jFixType"] = "ChangeSequence"
				else:
					json_val["jFixType"] = "Failed"
			else:
				json_val["jFixType"] = "Failed"

			# FINALISATION
			json_val["good"] = (json_val["vEnd"] != -1) and (json_val["jStart"] != -1)

			df.set_value(index, "cdr3fix." + gene_type, json.dumps(json_val))


	segments = pd.read_csv(segments_filepath, sep="\t")
	segments.columns = ["species", "gene", "segment", "id", "ref", "seq"]

	df = pd.read_csv(full_table, sep="\t")

	for index, row in df.iterrows():
		if (index+1) % 500 == 0:
			print(index + 1, "/", len(df))
		_fix_json(index, row, df, "alpha", segments)
		_fix_json(index, row, df, "beta", segments)

	df.to_csv(full_table, sep="\t", index=False)
	print("Done.")


if __name__ == "__main__":
	align_segments_and_write(sys.argv[1], sys.argv[2]) # vdjdb.full, segments.txt
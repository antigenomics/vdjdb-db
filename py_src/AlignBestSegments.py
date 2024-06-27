from __future__ import print_function
from __future__ import division

import sys
import json
import csv
import multiprocessing as mp

import pandas as pd

# Minimal number of nucleotides to make alignment
MIN_NUC_V = 5
MIN_NUC_J = 5
# Minimal difference among original segment and aligned segments using our aligner
MIN_DIFF_V = 2
MIN_DIFF_J = 2

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
    def _align_slice(codon_pos, gene_symbol, aminoacid, prev_aligned=[True] * 6):
        aligned = [False] * 6
        for i, codon in enumerate(CODONS[aminoacid]):
            aligned[i] = prev_aligned[i] & (codon[codon_pos] == gene_symbol)
        return sum(aligned) > 0, aligned

    score = 0

    flag = False

    for aa_pos, aa in enumerate(seq):
        if aa_pos * 3 < len(gene):
            flag, aligned = _align_slice(0, gene[aa_pos * 3], aa)

            if flag:
                score += 1

                if (aa_pos * 3 + 1) < len(gene):
                    flag, aligned = _align_slice(1, gene[aa_pos * 3 + 1], aa, aligned)

                    if flag:
                        score += 1

                        if (aa_pos * 3 + 2) < len(gene):
                            flag, aligned = _align_slice(2, gene[aa_pos * 3 + 2], aa, aligned)

                            score += flag
        if not flag:
            break

    return score


def align_nuc_to_aa_rev(seq, gene):
    def _align_slice(codon_pos, gene_symbol, aminoacid, prev_aligned=[True] * 6):
        aligned = [False] * 6
        for i, codon in enumerate(CODONS[aminoacid]):
            aligned[i] = prev_aligned[i] & (codon[codon_pos] == gene_symbol)
        return sum(aligned) > 0, aligned

    score = 0

    for aa_pos, aa in enumerate(reversed(seq)):
        if aa_pos * 3 < len(gene):
            flag, aligned = _align_slice(2, gene[-(aa_pos * 3 + 1)], aa)

            if flag:
                score += 1

                if (aa_pos * 3 + 1) < len(gene):
                    flag, aligned = _align_slice(1, gene[-(aa_pos * 3 + 2)], aa, aligned)

                    if flag:
                        score += 1

                        if (aa_pos * 3 + 2) < len(gene):
                            flag, aligned = _align_slice(0, gene[-(aa_pos * 3 + 3)], aa, aligned)

                            score += flag
        if not flag:
            break

    return score


def fix_json(row):
    cdr3 = row[0]
    res_v_id = "null"
    res_v_score = -1
    res_j_id = "null"
    res_j_score = -1
    if not (type(cdr3) is float):
        seg_gene_type = row[1]
        species = row[2]

        # VARIABLE
        for _, seg_row in segments[(segments.species == species) & (segments.gene == seg_gene_type) & (
                segments.segment == "Variable")].iterrows():
            cur_score = align_nuc_to_aa(cdr3, seg_row["seq"][seg_row["ref"] - 3:])
            if cur_score > res_v_score:
                if cur_score >= MIN_NUC_V:
                    res_v_score = cur_score
                    res_v_id = seg_row["id"]
        # elif cur_score == res_v_score:
        # 	res_v_id = res_v_id + "," + seg_row["id"]

        # JOINING
        for _, seg_row in segments[(segments.species == species) & (segments.gene == seg_gene_type) & (
                segments.segment == "Joining")].iterrows():
            cur_score = align_nuc_to_aa_rev(cdr3, seg_row["seq"][:seg_row["ref"] + 4])
            if cur_score > res_j_score:
                if cur_score >= MIN_NUC_J:
                    res_j_score = cur_score
                    res_j_id = seg_row["id"]
        # elif cur_score == res_j_score:
        # 	res_j_id = res_j_id + "," + seg_row["id"]

    return res_v_id, res_v_score // 3, res_j_id, res_j_score // 3


def update_segments(index, old_row, new_row, gene_type, single_col, stats, df):
    seg_gene_type = ""
    if gene_type == ".alpha":
        seg_gene_type = "TRA"
    elif gene_type == ".beta":
        seg_gene_type = "TRB"
    else:
        print("Error: unknown gene type", gene_type)
        return 0

    if single_col:
        gene_type = ""

    if not pd.isnull(row["cdr3fix" + gene_type]):
        json_val = json.loads(row["cdr3fix" + gene_type])

        if new_row[1] != -1:
            if new_row[0] == json_val["vId"]:
                if new_row[1] > json_val["vEnd"]:
                    json_val["oldVEnd"] = json_val["vEnd"]
                    json_val["vEnd"] = new_row[1]

                    json_val["oldVFixType"] = json_val["vFixType"]
                    json_val["vFixType"] = "Realign"
            elif (new_row[1] - json_val["vEnd"]) >= MIN_DIFF_V:
                json_val["oldVId"] = json_val["vId"]
                json_val["vId"] = new_row[0]

                json_val["oldVEnd"] = json_val["vEnd"]
                json_val["vEnd"] = new_row[1]

                json_val["oldVFixType"] = json_val["vFixType"]
                json_val["vFixType"] = "ChangeSegment"

        # if (new_row[1] - json_val["vEnd"]) >= MIN_DIFF_V:
        # 	json_val["oldVId"] = json_val["vId"]
        # 	json_val["vId"] = new_row[0]
        # 	json_val["oldVEnd"] = json_val["vEnd"]
        # 	json_val["vEnd"] = new_row[1]
        # 	json_val["oldVFixType"] = json_val["vFixType"]
        # 	if json_val["oldVId"] != json_val["vId"]:
        # 		temp_s = "v_" + str(json_val["oldVId"]) + "->" + json_val["vId"] + ":" + str(json_val["oldVEnd"]) + "->" + str(json_val["vEnd"])
        # 		stats[temp_s] = stats.get(temp_s, 0) + 1

        # 	if json_val["vId"] == json_val["oldVId"]:
        # 		json_val["vFixType"] = "Realign"
        # 	else:
        # 		json_val["vFixType"] = "ChangeSegment"

        if new_row[3] != -1:
            if new_row[2] == json_val["jId"]:
                if json_val["jStart"] > (len(row["cdr3" + gene_type]) - new_row[3]):
                    json_val["oldJStart"] = json_val["jStart"]
                    json_val["jStart"] = len(row["cdr3" + gene_type]) - new_row[3]

                    json_val["oldJFixType"] = json_val["jFixType"]
                    json_val["jFixType"] = "Realign"
            elif json_val["jStart"] - (len(row["cdr3" + gene_type]) - new_row[3]) >= MIN_DIFF_J:
                json_val["oldJId"] = json_val["jId"]
                json_val["jId"] = new_row[2]

                json_val["oldJStart"] = json_val["jStart"]
                json_val["jStart"] = len(row["cdr3" + gene_type]) - new_row[3]

                json_val["oldJFixType"] = json_val["jFixType"]
                json_val["jFixType"] = "ChangeSegment"

        # if (new_row[3] != -1) and ((json_val["jStart"] == -1) or (json_val["jStart"] - (len(row["cdr3" + gene_type]) - new_row[3]) >= MIN_DIFF_J)):
        # 	json_val["oldJId"] = json_val["jId"]
        # 	json_val["jId"] = new_row[2]
        # 	json_val["oldJStart"] = json_val["jStart"]
        # 	json_val["jStart"] = len(row["cdr3" + gene_type]) - new_row[3]
        # 	json_val["oldJFixType"] = json_val["jFixType"]
        # 	if json_val["oldJId"] != json_val["jId"]:
        # 		temp_s = "j_" + str(json_val["oldJId"]) + "->" + json_val["jId"] + ":" + str(json_val["oldJStart"]) + "->" + str(json_val["jStart"])
        # 		stats[temp_s] = stats.get(temp_s, 0) + 1

        # 	if json_val["jId"] == json_val["oldJId"]:
        # 		json_val["jFixType"] = "Realign"
        # 	else:
        # 		json_val["jFixType"] = "ChangeSegment"

        # FINALISATION
        json_val["good"] = (json_val["vEnd"] != -1) and (json_val["jStart"] != -1)

        df.set_value(index, "cdr3fix" + gene_type, json.dumps(json_val, sort_keys=True))

        if gene_type == "":
            gene_type = ".segm"
        df.set_value(index, "v" + gene_type, json_val["vId"])
        df.set_value(index, "j" + gene_type, json_val["jId"])

    return df


def json2tuples(df, cdr3seqcol, seg_gene_type):
    # json.loads(row["cdr3fix" + gene_type])
    return [(row[cdr3seqcol], seg_gene_type, row["species"]) for _, row in df.iterrows()]


def json2tuples2(df, cdr3seqcol):
    # json.loads(row["cdr3fix" + gene_type])
    return [(row[cdr3seqcol], row["gene"], row["species"]) for _, row in df.iterrows()]


def process_all(full_db, db):
    segments = pd.read_csv("../res/segments.txt", sep="\t")
    segments.columns = ["species", "gene", "segment", "id", "ref", "seq"]

    pool = mp.Pool(mp.cpu_count())

    print("-- processing the vdjdb_full.txt table")
    a_js = json2tuples(full_db, "cdr3.alpha", "TRA")
    b_js = json2tuples(full_db, "cdr3.beta", "TRB")
    a_js = pool.map(fix_json, a_js)
    b_js = pool.map(fix_json, b_js)

    print("-- writing")
    stats = {}
    for index, row in full_db.iterrows():
        full_db = update_segments(index, row, a_js[index], ".alpha", False, stats, df=full_db)
        full_db = update_segments(index, row, b_js[index], ".beta", False, stats, df=full_db)
    with open("../tmp/stats.full.txt", "w") as file:
        for s, v in sorted(stats.items(), key=lambda x: x[1], reverse=True):
            print(s, v, sep="\t", file=file)

    full_db.to_csv('../database/vdjdb.txt', sep="\t", index=False, quoting=csv.QUOTE_NONE)


	print("-- processing the vdjdb.txt table")
	ab_js = json2tuples2(db, "cdr3")
	ab_js = pool.map(fix_json, ab_js)

	print("-- writing")
	stats = {}
	for index, row in db.iterrows():
		db = update_segments(index, row, ab_js[index], ".alpha" if row["gene"] == "TRA" else ".beta", True, stats, df=db)
	with open("../tmp/stats.txt", "w") as file:
		for s, v in sorted(stats.items(), key=lambda x: x[1], reverse=True):
			print(s, v, sep="\t", file=file)

	db["method"] = db["method"].apply(json.loads).apply(json.dumps, sort_keys=True)
	db["meta"] = db["meta"].apply(json.loads).apply(json.dumps, sort_keys=True)
	db.to_csv('../database/vdjdb.txt', sep="\t", index=False, quoting=csv.QUOTE_NONE)


if __name__ == "__main__":
    # print(align_nuc_to_aa("CASSYLPGQGDHYSNQPQH", "TGTGCCAGCAGCTTAGG"))
    # # align_segments_and_write(sys.argv[1], sys.argv[2], sys.argv[3]) # vdjdb_full.txt, vdjdb.txt, segments.txt
    full_table = sys.argv[1]
    table = sys.argv[2]
    segments_filepath = sys.argv[3]

    segments = pd.read_csv(segments_filepath, sep="\t")
    segments.columns = ["species", "gene", "segment", "id", "ref", "seq"]

    pool = mp.Pool(mp.cpu_count())
    # pool = mp.Pool(1)

    print("-- processing the vdjdb_full.txt table")
    df = pd.read_csv(full_table, sep="\t")
    a_js = json2tuples(df, "cdr3.alpha", "TRA")
    b_js = json2tuples(df, "cdr3.beta", "TRB")
    a_js = pool.map(fix_json, a_js)
    b_js = pool.map(fix_json, b_js)

    print("-- writing")
    stats = {}
    for index, row in df.iterrows():
        update_segments(index, row, a_js[index], ".alpha", False, stats)
        update_segments(index, row, b_js[index], ".beta", False, stats)
    with open("../tmp/stats.full.txt", "w") as file:
        for s, v in sorted(stats.items(), key=lambda x: x[1], reverse=True):
            print(s, v, sep="\t", file=file)

    df.to_csv(full_table, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    print("-- processing the vdjdb.txt table")
    df = pd.read_csv(table, sep="\t")
    ab_js = json2tuples2(df, "cdr3")
    ab_js = pool.map(fix_json, ab_js)

    print("-- writing")
    stats = {}
    for index, row in df.iterrows():
        update_segments(index, row, ab_js[index], ".alpha" if row["gene"] == "TRA" else ".beta", True, stats)
    with open("../tmp/stats.txt", "w") as file:
        for s, v in sorted(stats.items(), key=lambda x: x[1], reverse=True):
            print(s, v, sep="\t", file=file)

    df["method"] = df["method"].apply(json.loads).apply(json.dumps, sort_keys=True)
    df["meta"] = df["meta"].apply(json.loads).apply(json.dumps, sort_keys=True)

    df.to_csv(table, sep="\t", index=False, quoting=csv.QUOTE_NONE)

def simplify_segment_name(segment_name: str) -> list:
    """
    simplifies segment name
    :param segment_name:
    :return: list of simplified version of segment name
    """
    no_allele = segment_name.split('*')[0]
    return [no_allele, segment_name.split('-')[0]]


translate_dict = {'TTT': 'F',
                  'TTC': 'F',
                  'TTA': 'L',
                  'TTG': 'L',
                  'TCT': 'S',
                  'TCC': 'S',
                  'TCA': 'S',
                  'TCG': 'S',
                  'TAT': 'Y',
                  'TAC': 'Y',
                  'TAA': '*',
                  'TAG': '*',
                  'TGT': 'C',
                  'TGC': 'C',
                  'TGA': '*',
                  'TGG': 'W',
                  'CTT': 'L',
                  'CTC': 'L',
                  'CTA': 'L',
                  'CTG': 'L',
                  'CCT': 'P',
                  'CCC': 'P',
                  'CCA': 'P',
                  'CCG': 'P',
                  'CAT': 'H',
                  'CAC': 'H',
                  'CAA': 'Q',
                  'CAG': 'Q',
                  'CGT': 'R',
                  'CGC': 'R',
                  'CGA': 'R',
                  'CGG': 'R',
                  'ATT': 'I',
                  'ATC': 'I',
                  'ATA': 'I',
                  'ATG': 'M',
                  'ACT': 'T',
                  'ACC': 'T',
                  'ACA': 'T',
                  'ACG': 'T',
                  'AAT': 'N',
                  'AAC': 'N',
                  'AAA': 'K',
                  'AAG': 'K',
                  'AGT': 'S',
                  'AGC': 'S',
                  'AGA': 'R',
                  'AGG': 'R',
                  'GTT': 'V',
                  'GTC': 'V',
                  'GTA': 'V',
                  'GTG': 'V',
                  'GCT': 'A',
                  'GCC': 'A',
                  'GCA': 'A',
                  'GCG': 'A',
                  'GAT': 'D',
                  'GAC': 'D',
                  'GAA': 'E',
                  'GAG': 'E',
                  'GGT': 'G',
                  'GGC': 'G',
                  'GGA': 'G',
                  'GGG': 'G'
                  }


def translate_linear(seq, reverse=False) -> str:
    """
    translates nucleotide sequence to amino acid sequence
    :param seq: nucleotide sequence to be translated
    :param reverse: if True reverse the sequence to be translated
    :return: translated sequence
    """
    aa_seq = ""
    if reverse:
        seq = seq[len(seq) % 3:]

    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if 'N' in codon and len(codon) == 3 and codon not in translate_dict.keys():
            aa_seq += 'X'
        elif codon not in translate_dict.keys() and len(codon) != 2:
            aa_seq += '?'
        else:
            aa_seq += translate_dict[codon]

    return aa_seq

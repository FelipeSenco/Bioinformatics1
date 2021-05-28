# DNA Toolkit:
from typing import Counter
from structures import *

# Check the sequence to make sure it is a DNA String:
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

# Count nucleotides in a DNA sequence:
def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

# Transcript DNa sequence to RNA:
def transcription(seq):
    """DNA -> RNA transcription. Replacing Thymine with Uracil"""
    return seq.replace("T", "U")

# Getting reverse complement DNA sequence:
def reverse_complement(seq):
    """Swapping adenine with thymine and guanine with cytosine. Reversing newly generated string"""
    # return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]
    # or
    # Pythonic approach. A little bit faster solution.
    mapping = str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)[::-1]

# Guanine and Cytocine content calculation:
def gc_content(seq):
    """GC Content in a DNA/RNA sequence"""
    return round((seq.count('C') + seq.count('G')) / len(seq) * 100)

def gc_content_subsec(seq, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. k=20 by default"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i: i + k]
        res.append(gc_content(subseq))
    return res

#translation function:
def translate_seq(seq, init_pos=0):
    """Translate  a DNA sequence into an aminoacid sequence"""
    return [DNA_Codons[seq[pos: pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]

#Condon usage function:
def condon_usage(seq, aminoacid):
    """Provides the frequency of each condon encoding a given aminoacid in a DNA sequence"""
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i: i + 3]] == aminoacid:
            tmpList.append(seq[i: i + 3])

    freqDict = dict(Counter(tmpList))
    totalWight = sum(freqDict.values())
    for codon in freqDict:
        freqDict[codon] = round(freqDict[codon] / totalWight, 2)
    return freqDict
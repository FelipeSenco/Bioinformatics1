#DNA Toolkit:
from structures import *

#Check the sequence to make sure it is a DNA String:
def validateSeq(dna_seq):
    tmpseq  = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

#Count nucleotides in a DNA sequence:
def countNucFrequency(seq):
    tmpFreqDict = {"A": 0,"C": 0,"G": 0,"T": 0}
    for nuc in seq:
       tmpFreqDict[nuc] += 1
    return tmpFreqDict

#Transcript DNa sequence to RNA:
def transcription(seq):
    """DNA -> RNA transcription. Replacing Thymine with Uracil"""
    return seq.replace("T", "U")

#Getting reverse complement DNA sequence:
def reverse_complement(seq):
    """Swapping adenine with thymine and guanine with cytosine. Reversing newly generated string"""
    return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]
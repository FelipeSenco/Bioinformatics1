#DNA Toolset/Code testing file
from DNAToolkit import *
import random
from utilities import colored

#Creating a random DNA sequence for test:
rndDNAStr = ''.join([random.choice(Nucleotides) 
                    for n in range(50)])


DNAStr = validateSeq(rndDNAStr)
RNAStr = transcription(DNAStr)

print(f'\nSequence: {colored(DNAStr)}\n')
print(f'[1] + Sequence length: {len(DNAStr)}\n')
print(f'[2] + Nucleotide frequency: {countNucFrequency(DNAStr)}\n')
print(f'[3] + DNA/RNA transcription: {colored(RNAStr)}\n')

print(f"[4] + DNA String + Receverse Complement:\n5' {(DNAStr)} 3'")
print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3' {(reverse_complement(DNAStr))} 5'\n")

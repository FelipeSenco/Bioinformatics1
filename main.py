#DNA Toolset/Code testing file
from DNAToolkit import *
import random
from utilities import colored

#Creating a random DNA sequence for test:
rndDNAStr = ''.join([random.choice(Nucleotides) 
                    for n in range(50)])


DNAStr = validateSeq(rndDNAStr)
RNAStr = transcription(DNAStr)

print(f'\nSequence: {(DNAStr)}\n')
print(f'[1] + Sequence length: {len(DNAStr)}\n')
print(f'[2] + Nucleotide frequency: {countNucFrequency(DNAStr)}\n')
print(f'[3] + DNA/RNA transcription: {(RNAStr)}\n')

print(f"[4] + DNA String + Complement + Receverse Complement:\n5' {(DNAStr)} 3'")
print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3' {(reverse_complement(DNAStr)[::-1])} 5' [Complement]")
print(f"5' {(reverse_complement(DNAStr))} 3' [Rev. Complement]\n")

print(f'[5] + GC Content: {gc_content(DNAStr)}%\n')
print(f'[6] + GC Content in Subsection k=5: {gc_content_subsec(DNAStr, k=5)}\n')

print(f'[7] + Aminoacids Sequence from DNA: {translate_seq(DNAStr, 0)}\n')
print(f"[8] + Condon frequency (L): {condon_usage(DNAStr, 'L')}")

print("[9] + Reading_frames:")
for frame in gen_reading_frames(DNAStr):
    print(frame)

print("\n[10] + All proteins in 6 open reading frames:")
for prot in all_proteins_from_orfs(DNAStr, startReadPos=0, endReadPos=0, ordered=True):
    print(f'{prot}')
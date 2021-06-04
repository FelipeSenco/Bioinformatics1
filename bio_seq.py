from bio_structures import DNA_Codons, DNA_Nucleotides
import random
from typing import Counter, Mapping, OrderedDict 

class bio_seq:
    """DNA sequence class. Default value: ATCG, DNA, No label"""

    def __init__(self, seq="ATCG", seq_type="DNA", label="No label"):
        """Sequence initialization, validation"""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"

    def __validate(self):
        """Check the sequence to make sure it is a valid DNA string"""
        return set(DNA_Nucleotides).issuperset(self.seq)

    def get_seq_biotype(self):
        """Return sequence type"""
        return self.seq_type

    def get_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}" 

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(DNA_Nucleotides)
                        for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    def nucleotide_frequency(self):
        """Returns a dictionary with each nucleotide count"""
        return dict(Counter(self.seq))
    
    def transcription(self):
        """DNA -> RNA Transcription"""
        return self.seq.replace("T", "U")

    def reverse_complement(self):
        """Returns the reverse complement of the DNA sequence"""
        mapping = str.maketrans("ATCG", "TAGC")
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        """GC Content in a DNA/RNA sequence"""
        return round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100)
    
    def gc_content_subsec(self, k=20):
        """GC Content in a DNA/RNA sub-sequence length k. k=20 by default"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i: i + k]
            res.append(round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100))
        return res
    
    def translate_seq(self, init_pos=0):
        """Translate  a DNA sequence into an aminoacid sequence"""
        return [DNA_Codons[self.seq[pos: pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
    
    def condon_usage(self, aminoacid):
        """Provides the frequency of each condon encoding a given aminoacid in a DNA sequence"""
        tmpList = []
        for i in range(0, len(self.seq) - 2, 3):
            if DNA_Codons[self.seq[i: i + 3]] == aminoacid:
                tmpList.append(self.seq[i: i + 3])

        freqDict = dict(Counter(tmpList))
        totalWight = sum(freqDict.values())
        for codon in freqDict:
            freqDict[codon] = round(freqDict[codon] / totalWight, 2)
        return freqDict

    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including the reverse complement"""
        frames = []
        frames.append(self.translate_seq(init_pos=0))
        frames.append(self.translate_seq(init_pos=1))
        frames.append(self.translate_seq(init_pos=2))

        tmp_reverse_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_reverse_seq.translate_seq(init_pos=0))
        frames.append(tmp_reverse_seq.translate_seq(init_pos=1))
        frames.append(tmp_reverse_seq.translate_seq(init_pos=2))
        del tmp_reverse_seq

        return frames

    def proteins_from_rf(self, aa_seq):
        """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                #STOP accumulating aminoacids if _ - STOP was found
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                # START accumulating amino acids if M - START was found
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""
        """Protein Seach DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
        """API can be used to pull protein info"""
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()
        
        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res
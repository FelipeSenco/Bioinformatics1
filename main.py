#DNA Toolset/Code testing file

from DNAToolkit import condon_usage
from bio_seq import bio_seq

test_dna = bio_seq()
test_dna.generate_rnd_seq(length=50)

print(test_dna.get_seq_info())
print(test_dna.nucleotide_frequency())
print(test_dna.transcription())
print(test_dna.reverse_complement())
print(test_dna.gc_content())
print(test_dna.gc_content_subsec())
print(test_dna.translate_seq())
print(test_dna.condon_usage('L'))
print(test_dna.gen_reading_frames())
print(test_dna.proteins_from_rf(['A', 'M', 'S', 'S', 'T', 'C', 'R', 'S', 'R', 'V', 'P', '_', 'S', 'W', 'T', 'N']))
print(test_dna.all_proteins_from_orfs())
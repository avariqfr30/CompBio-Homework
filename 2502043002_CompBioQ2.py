# Define a dictionary that maps DNA codons to mRNA codons
dna_to_mrna = {
    'A': 'U',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

# Define a dictionary that maps mRNA codons to amino acids (sideways layout)
mrna_to_amino_acid = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

#For number 1 :
def dna_to_mrna_translation(dna_sequence):
    mrna_sequence = ''.join([dna_to_mrna.get(base.upper(), 'Unknown') for base in dna_sequence])
    return mrna_sequence

def mrna_to_amino_acid_translation(mrna_sequence):
    amino_acid_sequence = []

    for i in range(0, len(mrna_sequence), 3):
        codon = mrna_sequence[i:i+3]
        if codon in mrna_to_amino_acid:
            amino_acid_sequence.append(mrna_to_amino_acid[codon])
        else:
            amino_acid_sequence.append('Unknown')
    return amino_acid_sequence

#For number 2 :
def count_rna_codon_frequency(dna_sequence, target_amino_acid=None):
    mrna_sequence = dna_to_mrna_translation(dna_sequence)
    amino_acid_sequence = mrna_to_amino_acid_translation(mrna_sequence)

    codon_frequency = {}

    for i, amino_acid in enumerate(amino_acid_sequence):
        if target_amino_acid is None or amino_acid == target_amino_acid:
            codon = mrna_sequence[i*3:i*3+3]
            if codon in codon_frequency:
                codon_frequency[codon] += 1
            else:
                codon_frequency[codon] = 1

    return codon_frequency

# (Number 1 related) Example DNA sequence
dna_sequence = "attacgattagc"

# (Number 2 related) Specify the target amino acid (or set it to None)
target_amino_acid = None # Change this to the desired amino acid or leave it as None

# (Number 2 related) Calculate the frequency of all RNA codons
frequency = count_rna_codon_frequency(dna_sequence, target_amino_acid)

# (Number 1 related) Print the DNA to mRNA translation
mrna_sequence = dna_to_mrna_translation(dna_sequence)
print("DNA to mRNA Translation:")
print("DNA Sequence:   " + dna_sequence.upper())  # Convert to uppercase for display
print("mRNA Sequence:  " + mrna_sequence)

# (Number 2 related) Print the amino acids of the DNA to mRNA sequence
amino_acids = mrna_to_amino_acid_translation(mrna_sequence)
print("\nAmino Acids:")
print(" - ".join(amino_acids))

if target_amino_acid is None:
    print("\nFrequencies of all RNA codons:")
    for codon, count in frequency.items():
        print(f"{codon}: {count}")
else:
    print(f"\nFrequencies of RNA codons encoding {target_amino_acid}:")
    for codon, count in frequency.items():
        print(f"{codon}: {count}")

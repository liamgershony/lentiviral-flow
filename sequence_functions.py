import os
import sys
import random

codons = {
       'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
       'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
       'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
       'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
       'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
       'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
       'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
       'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
       'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
       'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
       'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
       'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
       'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
       'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
       'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
       'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

# Generate dictionary of codons indexed by amino acid
aas = {}
for codon in codons:
    aa = codons[codon]
    if aa in aas:
        aas[aa].append(codon)
    else:
        aas[aa] = [codon]

def insert (source_str, insert_str, pos):
    return source_str[:pos]+insert_str+source_str[pos:]

def complement( seq ):
    seq = seq.upper()
    seq.replace('A','t')
    seq.replace('T','a')
    seq.replace('G','c')
    seq.replace('C','g')
    return seq.upper()

def translate(seq, frame = 0):
    seq = seq.upper()[frame:]
    protein =""
    for i in range(0, int(len(seq)/3)*3, 3): 
        codon = seq[i:i + 3]
        protein+= codons[codon]
    return protein

def find_subsequence_instances( sequence, subsequence ):
    start = 0
    seqlen = len(sequence)
    locs = []
    while subsequence in sequence[start:]:
        i = sequence.index(subsequence,start)
        locs.append(i)
        start = i+1
        if start >= seqlen:
            break
    return locs

def find_restriction_sites( restriction_sites, sequence ): 
    seqlen = len(sequence)
    locs = {}
    for s in restriction_sites:
        locs[s] = find_subsequence_instances( sequence, restriction_sites[s] )
    return locs

def read_gb( filename ):
    fin = open(filename,'r')
    record = False
    seq = ''
    for line in fin:
        if record:
            segments = line.lstrip().rstrip().split()[1:]
            for s in segments:
                seq += s
        elif 'ORIGIN' in line:
            record = True
        elif '//' == line[:2]:
            break
    fin.close()
    return seq
    #return translate(seq)

def check_codon_optimization( fseq, pseq ):

    (faa,paa) = get_intervening_sequence(fseq,pseq)

    print (faa)
    print (paa + '\n')

    faa = translate(faa)
    paa = translate(paa)

    print (faa)
    print (paa)

    if faa != paa:
        return False
    return True

# Randomly swap out one of the codons that intersect with
# the the taboo subsequence at the taboo index
def swap_codons( s, taboo, taboo_index ):
    # The index of the first nucleotide in the first codon that overlaps the taboo sequence
    first = 3*int(taboo_index/3)
    # The index of the first nucleotide in the last codon that overlaps the taboo sequence
    last = 3*int((taboo_index+len(taboo))/3)
    # The indices of the first nucleotides of the codons that overlap the taboo sequence
    affected_codon_indices = range(first, last+3, 3)
    # Randomly choose a codon that overlaps the taboo sequence
    codon_index = random.choice(affected_codon_indices)
    codon = s[codon_index:(codon_index+3)].upper()
    aa = codons[codon]
    # Cycle through possible codons that share the same amino acid and choose the next one
    for new_codon in sorted(aas[aa]):
        if new_codon != codon:
            break
    # Return the sequence with the offending codon switched
    return s[:codon_index] + new_codon + s[(codon_index+3):]

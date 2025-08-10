from Bio.Seq import Seq
from Bio.SeqUtils import GC, MeltingTemp as mt
import re

def nucleotide_count(seq):
    return { "A": seq.count("A"), "T": seq.count("T"), "G": seq.count("G"), "C": seq.count("C") }

def gc_content(seq):
    return GC(seq)

def complement(seq):
    return str(Seq(seq).complement())

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_repeats(seq, length=4):
    repeats = {}
    for i in range(len(seq) - length + 1):
        sub = seq[i:i+length]
        repeats[sub] = repeats.get(sub, 0) + 1
    return {k: v for k, v in repeats.items() if v > 1}

def find_orfs(seq):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    for i in range(len(seq)):
        if seq[i:i+3] == start_codon:
            for j in range(i+3, len(seq), 3):
                if seq[j:j+3] in stop_codons:
                    orfs.append(seq[i:j+3])
                    break
    return orfs

def melting_temp(seq):
    return mt.Tm_Wallace(seq)

def detect_cpg_islands(seq):
    islands = []
    for i in range(len(seq) - 200):
        window = seq[i:i+200]
        if GC(window) > 50 and window.count("CG") > 0.06 * len(window):
            islands.append((i, i+200))
    return islands

def primer_design(seq, length=20):
    primer = seq[:length]
    return {
        "primer": primer,
        "gc_content": GC(primer),
        "melting_temp": mt.Tm_Wallace(primer)
    }

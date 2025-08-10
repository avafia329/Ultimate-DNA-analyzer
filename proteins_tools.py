from Bio.Seq import Seq
from Bio.PDB import PDBParser, MMCIFParser
import io

def translate_dna(seq):
    return str(Seq(seq).translate(to_stop=True))

def simple_protein_structure(protein_seq):
    # Fake 3D coordinates for visualization demo
    coords = [{"x": i, "y": ord(aa) % 10, "z": (ord(aa) * 2) % 10} for i, aa in enumerate(protein_seq)]
    return coords

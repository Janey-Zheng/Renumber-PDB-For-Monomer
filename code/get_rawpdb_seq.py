# -*- coding: utf-8 -*-
import os
from pathlib import Path
from Bio.PDB import PDBParser

# Input folder containing .pdb files
STRUCT_DIR = r"F:\Data_0618\AlphaFold"

# Output folder for sequence text files
OUT_DIR = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\sequences\str_seq"
os.makedirs(OUT_DIR, exist_ok=True)

# Convert 3-letter to 1-letter amino acid code
def three_to_one(resname):
    mapping = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
        "GLU":"E","GLN":"Q","GLY":"G","HIS":"H","ILE":"I",
        "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
        "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
        "SEC":"U","PYL":"O"
    }
    return mapping[resname.upper()]

# Extract sequence and residue numbers per chain from a PDB file
def extract_chain_info(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = structure[0]  # first model

    chain_info = []
    for chain in model:
        residues = []
        res_ids = []
        for res in chain:
            hetfield, seqid, icode = res.id  # e.g. (' ', 42, 'A')
            # 跳过HETATM和水分子：标准氨基酸hetfield应为空格' '
            if hetfield.strip():  # 'H_'等；水通常为 'W'
                continue
            if "CA" not in res:  # 跳过不完整残基
                continue
            try:
                one_letter = three_to_one(res.resname)
            except KeyError:
                continue
            residues.append(one_letter)
            # 带上插入码（如果存在）
            res_ids.append(f"{seqid}{icode.strip()}" if icode.strip() else str(seqid))
        if residues:
            chain_info.append({
                "chain": chain.id,
                "sequence": "".join(residues),
                "ids": "-".join(res_ids)
            })
    return chain_info

# Write chain info into a text file
def write_txt(chain_info, outfile):
    with open(outfile, "w") as f:
        for item in chain_info:
            f.write(f"Chain: {item['chain']}\n")
            f.write(f"Sequence: {item['sequence']}\n")
            f.write(f"Auth Seq IDs: {item['ids']}\n\n")

def main():
    for root, dirs, files in os.walk(STRUCT_DIR):
        for fn in files:
            if fn.lower().endswith(".pdb"):
                pdb_file = os.path.join(root, fn)
                out_file = os.path.join(OUT_DIR, Path(fn).stem + ".txt")
                try:
                    chain_info = extract_chain_info(pdb_file)
                    write_txt(chain_info, out_file)
                    print(f"Processed: {pdb_file} -> {out_file}")
                except Exception as e:
                    print(f"Failed on {pdb_file}: {e}")

if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
import os
from pathlib import Path
from tqdm import tqdm
from Bio.PDB import PDBParser

# Input root: each subfolder contains *_pocket_clean.pdb and *_protein_clean_singlechain.pdb
IN_ROOT = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\PDBbind\PDBbind_structure"

# Output root: one folder per structure (same name as subfolder under IN_ROOT)
OUT_ROOT = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\sequences\pdbbind_seq"
os.makedirs(OUT_ROOT, exist_ok=True)


# Map 3-letter residue names to 1-letter codes (standard + SEC/ PYL)
AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLU":"E","GLN":"Q","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "SEC":"U","PYL":"O"
}


# Extract per-chain sequence and residue numbers (auth_seq_id ~ seqid) from a PDB file
def extract_chain_info_from_pdb(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    model = structure[0]  # first model

    chain_blocks = []
    for chain in model:
        seq_chars = []
        id_list = []
        for res in chain:
            hetflag, seqid, icode = res.id  # (' ', 123, ' ')
            if hetflag.strip():     # skip HETATM residues
                continue
            if "CA" not in res:     # skip very incomplete residues
                continue
            resname = res.resname.upper()
            one = AA3_TO_1.get(resname)
            if not one:             # skip non-standard residues we can't map
                continue
            seq_chars.append(one)
            id_list.append(str(seqid))  # use seqid as auth_seq_id

        if seq_chars:
            chain_blocks.append({
                "chain": chain.id,
                "sequence": "".join(seq_chars),
                "ids": "-".join(id_list)
            })
    return chain_blocks


# Write a list of chain blocks to a .txt file in the requested format
def write_chain_txt(blocks, out_txt):
    with open(out_txt, "w") as f:
        for blk in blocks:
            f.write(f"Chain: {blk['chain']}\n")
            f.write(f"Sequence: {blk['sequence']}\n")
            f.write(f"Auth Seq IDs: {blk['ids']}\n\n")


# For one structure folder (e.g., 1a42), process pocket/protein files and write to OUT_ROOT/<name>/
def process_structure_folder(struct_dir):
    name = Path(struct_dir).name
    pocket = None
    protein = None
    # find files
    for fn in os.listdir(struct_dir):
        if fn.endswith("_pocket_clean_singlechain.pdb"):
            pocket = os.path.join(struct_dir, fn)
        elif fn.endswith("_protein_clean_singlechain.pdb"):
            protein = os.path.join(struct_dir, fn)

    # skip if neither file exists
    if not pocket and not protein:
        return False, f"No target files in {struct_dir}"

    # make output subfolder
    out_dir = os.path.join(OUT_ROOT, name)
    os.makedirs(out_dir, exist_ok=True)

    # process pocket
    if pocket:
        blocks = extract_chain_info_from_pdb(pocket)
        out_txt = os.path.join(out_dir, Path(pocket).stem + ".txt")
        write_chain_txt(blocks, out_txt)

    # process protein (single chain expected, but code supports any)
    if protein:
        blocks = extract_chain_info_from_pdb(protein)
        out_txt = os.path.join(out_dir, Path(protein).stem + ".txt")
        write_chain_txt(blocks, out_txt)

    return True, f"OK {name}"


# Walk through all subfolders under IN_ROOT and process them with a progress bar
def main():
    subdirs = [os.path.join(IN_ROOT, d) for d in os.listdir(IN_ROOT)
               if os.path.isdir(os.path.join(IN_ROOT, d))]
    for d in tqdm(subdirs, desc="Generating PDBbind sequences", unit="dir"):
        try:
            ok, msg = process_structure_folder(d)
        except Exception as e:
            ok, msg = False, f"Error {d}: {e}"
        if not ok:
            print(msg)


if __name__ == "__main__":
    main()

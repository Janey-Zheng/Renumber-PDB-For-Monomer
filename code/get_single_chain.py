# -*- coding: utf-8 -*-
import os
from pathlib import Path

# Root folder of PDBbind_structure
ROOT_DIR = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\PDBbind\PDBbind_structure"


# Remove water molecules (HOH) and hetero atoms (HETATM) from a PDB file
def clean_pdb(infile, outfile):
    with open(infile, "r") as fin, open(outfile, "w") as fout:
        for line in fin:
            if line.startswith("HETATM"):
                continue
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                if resname == "HOH":
                    continue
                fout.write(line)
            elif line.startswith("TER") or line.startswith("END"):
                fout.write(line)


# Keep only the first chain in a cleaned PDB file (robust to short lines / blank chain IDs)
# Keep only the alphabetically first chain in a cleaned PDB file
def keep_first_chain(infile, outfile):
    # Step 1: collect all chain IDs
    chain_ids = set()
    with open(infile, "r") as fin:
        for line in fin:
            if line.startswith("ATOM"):
                if len(line) > 21:
                    chain_ids.add(line[21].strip() or " ")
    if not chain_ids:
        return  # nothing to keep

    # Step 2: pick the alphabetically first chain
    first_chain = sorted(chain_ids)[0]

    # Step 3: write only this chain
    with open(infile, "r") as fin, open(outfile, "w") as fout:
        for line in fin:
            rec = line[:6]
            if rec.startswith("ATOM") or rec.startswith("TER"):
                chain_id = line[21] if len(line) > 21 else " "
                if chain_id.strip() == first_chain:
                    fout.write(line)
            elif rec.startswith("END"):
                fout.write(line)
            else:
                # keep headers/remarks if present
                fout.write(line)


# Process one pdb file: clean -> keep first chain -> report paths
def process_one_pdb(pdb_path):
    cleaned = os.path.join(os.path.dirname(pdb_path), Path(pdb_path).stem + "_clean.pdb")
    single  = os.path.join(os.path.dirname(pdb_path), Path(pdb_path).stem + "_clean_singlechain.pdb")
    clean_pdb(pdb_path, cleaned)
    print(f"Cleaned: {pdb_path} -> {cleaned}")
    keep_first_chain(cleaned, single)
    print(f"Single chain: {cleaned} -> {single}")


# Traverse subfolders, process both *_pocket.pdb and *_protein.pdb with the same pipeline
def main():
    for root, _, files in os.walk(ROOT_DIR):
        for fn in files:
            if fn.endswith("_pocket.pdb") or fn.endswith("_protein.pdb"):
                infile = os.path.join(root, fn)
                process_one_pdb(infile)


if __name__ == "__main__":
    main()

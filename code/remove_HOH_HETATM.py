# -*- coding: utf-8 -*-
import os
from pathlib import Path

# Root folder of PDBbind_structure
ROOT_DIR = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\PDBbind\PDBbind_structure"


# Function: remove water molecules (HOH) and hetero atoms (HETATM) from a PDB file
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


# Function: traverse all subfolders and clean *_pocket.pdb and *_protein.pdb files
def main():
    for root, dirs, files in os.walk(ROOT_DIR):
        for fn in files:
            if fn.endswith("_pocket.pdb") or fn.endswith("_protein.pdb"):
                infile = os.path.join(root, fn)
                outfile = os.path.join(root, Path(fn).stem + "_clean.pdb")
                clean_pdb(infile, outfile)
                print(f"Cleaned: {infile} -> {outfile}")


if __name__ == "__main__":
    main()

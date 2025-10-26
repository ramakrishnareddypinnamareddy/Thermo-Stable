#!/usr/bin/env python3
# coding: utf-8

import os
import sys

def repair(p):
    """Repair a PDB file using FoldX"""
    cmd = f"./foldx --command=RepairPDB --pdb={p}.pdb"
    os.system(cmd)

def free_energy_membrane(p):
    """Calculate free energy using FoldX Stability command"""
    cmd = f"./foldx --command=Stability --pdb={p}.pdb --output-dir=./foldx_res"
    os.system(cmd)  # run FoldX

    out_file = f'./foldx_res/{p}_0_ST.fxout'
    if not os.path.exists(out_file):
        raise FileNotFoundError(f"FoldX output file not found: {out_file}")

    with open(out_file, 'r') as f:
        lines = f.readlines()

    # Skip header, get first data line
    for line in lines:
        if line.startswith("Total"):
            continue
        parts = line.strip().split('\t')
        # Columns we need: BackHbond (1), SideHbond (2), Electro (5), Energy_SolvP (6)
        back_hbond = float(parts[1])
        side_hbond = float(parts[2])
        electro = float(parts[5])
        solv_p = float(parts[6])
        return back_hbond, side_hbond, electro, solv_p

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 foldx.py <protein1.pdb> <protein2.pdb>")
        sys.exit(1)

    pdb1, pdb2 = sys.argv[1], sys.argv[2]
    pdbs = [pdb1, pdb2]

    results = []

    for pdb_file in pdbs:
        protein_name = os.path.splitext(pdb_file)[0]
        print(f"Processing {protein_name}...")
        features = free_energy_membrane(protein_name)
        results.append(features)
        print(f"{protein_name} Features: {features}")

    # Output for run_classifier parsing
    print("\t".join(str(x) for x in results[0]) + "\t" + "\t".join(str(x) for x in results[1]))

if __name__ == "__main__":
    main()


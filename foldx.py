#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd

def repair(p):
    """Repair a PDB file using FoldX"""
    cmd = f"./foldx --command=RepairPDB --pdb={p}.pdb"
    os.system(cmd)  # just run the command, no return needed

def free_energy_membrane(p):
    """Calculate free energy using FoldX Stability command"""
    #p=p.split('/)[-1]
    #print(p)
    cmd = f"./foldx --command=Stability --pdb={p}.pdb --output-dir=./foldx_res"
    #print(cmd)
    os.system(cmd)  # run FoldX
    out_file = f'./foldx_res/{p}_0_ST.fxout'
    if not os.path.exists(out_file):
        raise FileNotFoundError(f"FoldX output file not found: {out_file}")
    with open(out_file, 'r') as f:
        lines = f.read().strip()
    values = lines.split('\t')
    return values


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 run_foldx.py <protein1.pdb> <protein2.pdb>")
        sys.exit(1)

    pdb1, pdb2 = sys.argv[1], sys.argv[2]
    pdbs = [pdb1, pdb2]

    df_columns = [
        'entry','Total','BackHbond','SideHbond','Energy_VdW','Electro','Energy_SolvP','Energy_SolvH',
        'Energy_vdwclash','Entropy_sidec','Entropy_mainc','water bonds','dummy','cis_bond','energy_torsion',
        'backbone_vdwclash','helix dipole','loop_entropy','disulfide','kn electrostatic','partial covalent interactions',
        'Energy_Ionisation','Entropy Complex','count'
    ]

    results = {}

    for pdb_file in pdbs:
        protein_name = os.path.splitext(pdb_file)[0]
        print(f"Processing {protein_name}...")

        # Run FoldX stability
        energy_values = free_energy_membrane(protein_name)
        results[protein_name] = energy_values

        # Only get Electro value (6th column, index 5)
        electro_value = energy_values[5]
        print(f"{protein_name} Electro: {electro_value}")

    return results

if __name__ == "__main__":
    main()



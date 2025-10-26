#!/usr/bin/env python
# coding: utf-8

from Bio.PDB import *
import sys
import pandas as pd
import os

# --- Amino acid mapping ---
AA = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
      'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

# --- Read CA residues from PDB ---
def PDB_Read(pdb_file):   
    with open(f'./{pdb_file}', 'r') as f:
        lines = f.readlines()
    ca_lines = [i for i in lines if len(i.split()) > 3 and i.split()[0] == 'ATOM' and i.split()[2] == 'CA']
    resnames = [split_pdb(i)[3] for i in ca_lines]
    res1 = [AA[i] for i in resnames if i in AA.keys()]
    return res1

# --- Split PDB line ---
def split_pdb(line):
    splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
    return splitted_line

# --- Round helper ---
def roundd(a):
    return round(float(a), 3)

# --- Main execution ---
def main():
    if len(sys.argv) != 3:
        print("Usage: python3 aaindex_script.py <protein1.pdb> <protein2.pdb>")
        sys.exit(1)

    pdb1, pdb2 = sys.argv[1], sys.argv[2]
    pdbs = [pdb1, pdb2]

    # Read AAIndex CSV (assumes first row is header, first column is residue)
    df = pd.read_csv('aaindex.csv')

    # Only two feature columns assumed
    feature_columns = df.columns[2:7]  # columns 2 and 3
    feature_names = [df.iloc[0, c] for c in range(2,7)]

    results = {}

    for pdb_file in pdbs:
        protein_name = os.path.splitext(os.path.basename(pdb_file))[0]
        residues = PDB_Read(pdb_file)

        feature_values = []
        for col in range(2, 7):
            keys = list(df.iloc[2:, 1])
            values = list(df.iloc[2:, col])
            aa_dict = {keys[i]: float(values[i]) for i in range(len(keys))}
            s = sum([aa_dict[res] for res in residues])
            avg = s / len(residues) if residues else 0
            feature_values.append(avg)

        results[protein_name] = feature_values

        # Print only the two feature values
        print(f"\n=== {protein_name} AAIndex Features ===")
        for fname, val in zip(feature_names, feature_values):
            print(f"{fname}: {roundd(val)}")

    return results

if __name__ == "__main__":
    main()


#!/usr/bin/env python
import sys
import subprocess
import pandas as pd
import re

# Define residue groups
CHARGED = {"ARG", "LYS", "ASP", "GLU", "HIS"}
APOLAR = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "PRO", "TRP"}
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
POLAR_UNCHARGED = {"SER", "THR", "ASN", "GLN", "CYS", "GLY"}

def run_script(script, pdb1, pdb2):
    result = subprocess.run(
        ["python3", script, pdb1, pdb2],
        capture_output=True,
        text=True
    )
    return result.stdout

def parse_sride_output(output):
    counts = re.findall(r"Total stabilizing residues: (\d+)", output)
    return [int(c) for c in counts]

def parse_foldx_output(output):
    electro_values = re.findall(r"Electro:\s*([-\d\.]+)", output)
    return [float(e) for e in electro_values]

def parse_aaindex_output(output):
    features = []
    for line in output.splitlines():
        if ":" in line:
            val = float(line.split(":")[1].strip())
            features.append(val)
    return [features[:2], features[2:]]

def parse_pdb_composition(pdb_file):
    """Count residue types from PDB file"""
    res_counts = {"length": 0, "charged": 0, "apolar": 0, "aromatic": 0, "polar_uncharged": 0}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and line[13:15].strip() == "CA":
                res = line[17:20].strip()
                res_counts["length"] += 1
                if res in CHARGED:
                    res_counts["charged"] += 1
                if res in APOLAR:
                    res_counts["apolar"] += 1
                if res in AROMATIC:
                    res_counts["aromatic"] += 1
                if res in POLAR_UNCHARGED:
                    res_counts["polar_uncharged"] += 1
    return res_counts

def classify(df):
    score_P1 = 0
    score_P2 = 0

    # Original features
    features = [("StabRes", True),
                ("Electro", False),
                ("Feature2", True),
                ("Feature3", True)]

    # New features
    features += [("Length", True),
                 ("Charged", True),
                 ("Apolar", True),
                 ("Aromatic", True),
                 ("PolarUncharged", False)]  # less polar uncharged is better

    for feature, higher_is_better in features:
        p1_val = df[f"P1_{feature}"].iloc[0]
        p2_val = df[f"P2_{feature}"].iloc[0]

        diff = p2_val - p1_val
        if diff == 0:
            continue
        if higher_is_better:
            if diff > 0:
                score_P2 += 1
            else:
                score_P1 += 1
        else:
            if diff < 0:
                score_P2 += 1
            else:
                score_P1 += 1

    # Decide thermostable protein
    if score_P1 > score_P2:
        print(f"{df['Protein1'].iloc[0]} is more thermostable")
    elif score_P2 > score_P1:
        print(f"{df['Protein2'].iloc[0]} is more thermostable")
    else:
        print("Unable to decide thermostable protein: tied scores")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 run_classifier.py <protein1.pdb> <protein2.pdb>")
        sys.exit(1)

    pdb1, pdb2 = sys.argv[1], sys.argv[2]

    # Run scripts
    stab_counts = parse_sride_output(run_script("SRIDE.py", pdb1, pdb2))
    electro_values = parse_foldx_output(run_script("foldx.py", pdb1, pdb2))
    aa_features = parse_aaindex_output(run_script("aaindex.py", pdb1, pdb2))

    # Compute composition features
    comp1 = parse_pdb_composition(pdb1)
    comp2 = parse_pdb_composition(pdb2)

    data = {
        "Protein1": [pdb1],
        "Protein2": [pdb2],
        "P1_StabRes": [stab_counts[0]],
        "P1_Electro": [electro_values[0]],
        "P1_Feature2": [aa_features[0][0]],
        "P1_Feature3": [aa_features[0][1]],
        "P1_Length": [comp1["length"]],
        "P1_Charged": [comp1["charged"]],
        "P1_Apolar": [comp1["apolar"]],
        "P1_Aromatic": [comp1["aromatic"]],
        "P1_PolarUncharged": [comp1["polar_uncharged"]],
        "P2_StabRes": [stab_counts[1]],
        "P2_Electro": [electro_values[1]],
        "P2_Feature2": [aa_features[1][0]],
        "P2_Feature3": [aa_features[1][1]],
        "P2_Length": [comp2["length"]],
        "P2_Charged": [comp2["charged"]],
        "P2_Apolar": [comp2["apolar"]],
        "P2_Aromatic": [comp2["aromatic"]],
        "P2_PolarUncharged": [comp2["polar_uncharged"]],
    }

    df = pd.DataFrame(data)
    classify(df)


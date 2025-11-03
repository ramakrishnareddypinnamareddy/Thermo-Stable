#!/usr/bin/env python3
import sys
import subprocess
import pandas as pd
import re
from Bio import SeqIO


import zipfile

with zipfile.ZipFile("protein_thermo_meso_classifier.zip", "r") as zip_ref:
    zip_ref.extractall(".")  # Extracts all contents to the current directory

import tarfile

# Path to your .tar.gz file
tar_path = "foldx.tar.xz"
import os

folder_name = "foldx_res"
os.makedirs(folder_name, exist_ok=True)  # Creates the folder if it doesn’t exist

# Extract all contents to the current directory
with tarfile.open(tar_path, "r:xz") as tar_ref:
    tar_ref.extractall(".")  # You can replace "." with a folder path if desired

    
# === CONFIG ===
CSV_FILE = "meso_thermo_sequences.csv"  # update to your CSV file name
MESO_COL = "Meso_Sequence"
THERMO_COL = "Thermo_Sequence"
# Manual three-letter to one-letter mapping
three_to_one_dict = {
    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
    'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
    'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y',
    'SEC':'U', 'PYL':'O', 'ASX':'B', 'GLX':'Z', 'XLE':'J', 'XAA':'X'
}

import pandas as pd

def save_protein_features(data, output_file="protein_features.csv"):
    """
    Convert a dictionary of protein features into a DataFrame with proper column names
    and save it as a CSV file.

    Parameters:
    - data (dict): Dictionary containing protein features.
    - output_file (str): Path to save the CSV file. Default is 'protein_features.csv'.

    Returns:
    - df (pd.DataFrame): The resulting DataFrame.
    """

    # Column mapping
    column_mapping = {
        "Protein1": "P1",
        "Protein2": "P2",
        "P1_Feature2": "P1_Apparent partial specific volume (Bull-Breese, 1974)",
        "P1_Feature3": "P1_Localized electrical effect (Fauchere et al., 1988)",
        "P1_Feature4": "P1_Short and medium range non-bonded energy per residue (Oobatake-Ooi, 1977)",
        "P1_Feature5": "P1_Total weighted degree of the graph (obtained by adding all the weights of",
        "P1_Feature6": "P1_pK-N (Fasman, 1976)",
        "P1_StabRes": "P1_Total_Stabilizing",
        "P1_Degree": "P1_Degree",
        "P1_Clustering": "P1_Clustering Coeff",
        "P1_ANN": "P1_ANN Degree",
        "P1_BackHbond": "P1_BackHbond",
        "P1_SideHbond": "P1_SideHbond",
        "P1_Electro": "P1_Electro",
        "P1_Energy_SolvP": "P1_Energy_SolvP",
        "P2_Feature2": "P2_Apparent partial specific volume (Bull-Breese, 1974)",
        "P2_Feature3": "P2_Localized electrical effect (Fauchere et al., 1988)",
        "P2_Feature4": "P2_Short and medium range non-bonded energy per residue (Oobatake-Ooi, 1977)",
        "P2_Feature5": "P2_Total weighted degree of the graph (obtained by adding all the weights of",
        "P2_Feature6": "P2_pK-N (Fasman, 1976)",
        "P2_StabRes": "P2_Total_Stabilizing",
        "P2_Degree": "P2_Degree",
        "P2_Clustering": "P2_Clustering Coeff",
        "P2_ANN": "P2_ANN Degree",
        "P2_BackHbond": "P2_BackHbond",
        "P2_SideHbond": "P2_SideHbond",
        "P2_Electro": "P2_Electro",
        "P2_Energy_SolvP": "P2_Energy_SolvP",
    }

    # Convert to DataFrame
    df = pd.DataFrame(data)

    # Rename columns
    df = df.rename(columns=column_mapping)

    # Save to CSV
    #df.to_csv(output_file, index=False)
    #print(f"DataFrame saved as {output_file}")

    return df


def three_to_one(resname):
    """
    Convert 3-letter residue code to 1-letter.
    Unknown residues return 'X'.
    """
    return three_to_one_dict.get(resname.upper(), 'X')

def get_fasta_sequence_from_pdb(pdb_file):
    """
    Extract protein sequence from ATOM lines of PDB.
    Returns a single concatenated sequence of all chains.
    """
    seq = []
    seen_residues = set()
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                chain_id = line[21].strip()
                resseq = int(line[22:26].strip())
                res_id = (chain_id, resseq)
                if res_id not in seen_residues:
                    seen_residues.add(res_id)
                    seq.append(three_to_one(resname))
    return "".join(seq)

def run_script(script, pdb1, pdb2):
    """Run an external Python script and return stdout."""
    result = subprocess.run(["python3", script, pdb1, pdb2],
                            capture_output=True, text=True)
    return result.stdout

def parse_sride_output(output):
    #print("FoldX")
    #print(output)
    """Extract stabilizing residues count."""
    counts = re.findall(r"Total stabilizing residues: (\d+)", output)
    return [int(c) for c in counts]

def parse_foldx_output(output):
    """
    Extract BackHbond, SideHbond, Electro, Energy_SolvP values.
    Output format from foldx.py: last line contains 8 numbers (tab-separated).
    """
    lines = output.strip().splitlines()
    # Take the last non-empty line
    for line in reversed(lines):
        line = line.strip()
        if not line:
            continue
        # check if it looks like feature numbers (contains tabs or multiple numbers)
        parts = line.replace(',', '.').split()  # replace commas if any
        if len(parts) == 8:
            p1_features = [float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])]
            p2_features = [float(parts[4]), float(parts[5]), float(parts[6]), float(parts[7])]
            return [p1_features, p2_features]

    raise ValueError(f"Unexpected FoldX output: {output}")


def parse_aaindex_output(output):
    """
    Extract AAIndex features for two proteins from output.
    Returns a list of two lists: [[protein1_features], [protein2_features]]
    """
    features = []
    current_protein_features = []
    
    for line in output.splitlines():
        line = line.strip()
        # Skip empty lines
        if not line:
            continue
        # Detect protein header
        if "AAIndex Features" in line:
            if current_protein_features:
                features.append(current_protein_features)
                current_protein_features = []
        # Extract numeric value after colon
        elif ":" in line:
            try:
                val = float(line.split(":")[1].strip())
                current_protein_features.append(val)
            except ValueError:
                # Some lines might not have valid float
                continue
    
    # Append last protein
    if current_protein_features:
        features.append(current_protein_features)
    
    return features

def parse_network_output(output):
    """
    Extract Degree, Clustering Coefficient, ANN Degree for both proteins.
    Returns: [[P1_degree, P1_clust, P1_ann], [P2_degree, P2_clust, P2_ann]]
    """
    vals = [float(x) for x in output.strip().split()]
    if len(vals) != 6:
        raise ValueError(f"Unexpected network output: {output}")
    return [vals[:3], vals[3:]]

def check_sequence_in_csv(seq1, seq2, csv_file):
    """Check if sequences exist in Meso_Sequence or Thermo_Sequence columns."""
    df = pd.read_csv(csv_file)
    seq1_match, seq2_match = None, None

    for seq, label in [(seq1, "Protein1"), (seq2, "Protein2")]:
        if seq in df[THERMO_COL].values:
            match = "Thermo"
        elif seq in df[MESO_COL].values:
            match = "Meso"
        else:
            match = None
        if label == "Protein1":
            seq1_match = match
        else:
            seq2_match = match

    return seq1_match, seq2_match

def classify(df):
    """Feature-based classification."""
    score_P1, score_P2 = 0, 0
    features = [
        ("StabRes", True),
        ("BackHbond", True),
        ("SideHbond", False),
        ("Electro", False),
        ("Energy_SolvP", True),
        ("Feature2", True),
        ("Feature3", True),
                ("Feature4", True),
        ("Feature5", True), ("Feature6", True), ("Degree", True),
        ("Clustering", True),
        ("ANN", True)
    ]
    for feature, higher_is_better in features:
        p1_val = df[f"P1_{feature}"].iloc[0]
        p2_val = df[f"P2_{feature}"].iloc[0]

        if p1_val == p2_val:
            continue

        if higher_is_better:
            if p2_val > p1_val:
                score_P2 += 1
            else:
                score_P1 += 1
        else:  # Electro: less negative is better
            if p2_val < p1_val:
                score_P2 += 1
            else:
                score_P1 += 1

    if score_P1 > score_P2:
        print(f"{df['Protein1'].iloc[0]} is more thermostable")
    elif score_P2 > score_P1:
        print(f"{df['Protein2'].iloc[0]} is more thermostable")
    else:
        p1_stab, p2_stab = df["P1_StabRes"].iloc[0], df["P2_StabRes"].iloc[0]
        if p1_stab > p2_stab:
            print(f"{df['Protein1'].iloc[0]} is more thermostable")
        elif p2_stab > p1_stab:
            print(f"{df['Protein2'].iloc[0]} is more thermostable")
        else:
            p1_elec, p2_elec = df["P1_Electro"].iloc[0], df["P2_Electro"].iloc[0]
            if p1_elec > p2_elec:
                print(f"{df['Protein1'].iloc[0]} is more thermostable")
            elif p2_elec > p1_elec:
                print(f"{df['Protein2'].iloc[0]} is more thermostable")
            else:
                print("Unable to decide thermostable protein: perfectly tied")
import pandas as pd
import pickle

def predict_protein_pairs(df, model_file):
    """
    Predict Mesophile/Thermophile (or stability) classes for P1 and P2 proteins,
    and print the results as a table: Protein | Prob_Thermo | Class
    """
    # === Load model ===
    with open(model_file, "rb") as f:
        model = pickle.load(f)

    # === Base features ===
    features = [
        "Apparent partial specific volume (Bull-Breese, 1974)",
        "Localized electrical effect (Fauchere et al., 1988)",
        "Short and medium range non-bonded energy per residue (Oobatake-Ooi, 1977)",
        "Total weighted degree of the graph (obtained by adding all the weights of",
        "pK-N (Fasman, 1976)",
        "Total_Stabilizing",
        "Degree",
        "Clustering Coeff",
        "BackHbond",
        "SideHbond",
        "Electro",
        "Energy_SolvP"
    ]

    # === Build feature matrices ===
    X_p1 = df[[f"P1_{f}" for f in features]].copy()
    X_p1.columns = features
    X_p2 = df[[f"P2_{f}" for f in features]].copy()
    X_p2.columns = features

    # === Predict ===
    p1_pred = model.predict(X_p1)
    p2_pred = model.predict(X_p2)

    # Probabilities (if available)
    if hasattr(model, "predict_proba"):
        p1_prob = model.predict_proba(X_p1)[:, 1]
        p2_prob = model.predict_proba(X_p2)[:, 1]
    else:
        p1_prob = [None] * len(df)
        p2_prob = [None] * len(df)

    # Map predictions to class labels
    p1_class = [ "Less Stable" if x==0 else "More Stable" for x in p1_pred ]
    p2_class = [ "Less Stable" if x==0 else "More Stable" for x in p2_pred ]

    # === Build long-format table for printing ===
    table = pd.DataFrame({
        "Protein": list(df["P1"]) + list(df["P2"]),
        "Prob_Thermo": list(p1_prob) + list(p2_prob),
        "Class": list(p1_class) + list(p2_class)
    })

    # Print table
    print("\n=== Prediction Table ===")
    print(table.to_string(index=False))

    return table

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 run_classifier.py <protein1.pdb> <protein2.pdb>")
        sys.exit(1)

    pdb1, pdb2 = sys.argv[1], sys.argv[2]

    # Extract sequences
    seq1 = get_fasta_sequence_from_pdb(pdb1)
    seq2 = get_fasta_sequence_from_pdb(pdb2)

    if not seq1 or not seq2:
        print("Error: Could not extract sequences from one or both PDB files.")
        sys.exit(1)

    # === Skip checking CSV and directly calculate ===
    print("Proceeding with feature-based classification...")
    print()

    print("Calculating Stabilizing residues...\n")
    sride_out = run_script("SRIDE.py", pdb1, pdb2)
    stab_counts = parse_sride_output(sride_out)

    print("Calculating Sequence-Based Features...\n")
    aaindex_out = run_script("aaindex.py", pdb1, pdb2)
    aa_features = parse_aaindex_output(aaindex_out)

    print("Calculating FoldX Features...\n")
    foldx_out = run_script("foldx.py", pdb1, pdb2)
    foldx_features = parse_foldx_output(foldx_out)

    print("Calculating Network Features...\n")
    network_out = run_script("network_features.py", pdb1, pdb2)
    network_features = parse_network_output(network_out)

    # Combine all features into DataFrame
    data = {
        "Protein1": [pdb1],
        "Protein2": [pdb2],
        "P1_StabRes": [stab_counts[0]],
        "P1_BackHbond": [foldx_features[0][0]],
        "P1_SideHbond": [foldx_features[0][1]],
        "P1_Electro": [foldx_features[0][2]],
        "P1_Energy_SolvP": [foldx_features[0][3]],
        "P2_StabRes": [stab_counts[1]],
        "P2_BackHbond": [foldx_features[1][0]],
        "P2_SideHbond": [foldx_features[1][1]],
        "P2_Electro": [foldx_features[1][2]],
        "P2_Energy_SolvP": [foldx_features[1][3]],
        "P1_Feature2": [aa_features[0][0]],
        "P1_Feature3": [aa_features[0][1]],
        "P2_Feature2": [aa_features[1][0]],
        "P2_Feature3": [aa_features[1][1]],
        "P1_Feature4": [aa_features[0][2]],
        "P1_Feature5": [aa_features[0][3]],
        "P2_Feature4": [aa_features[1][2]],
        "P2_Feature5": [aa_features[1][3]],
        "P1_Feature6": [aa_features[0][4]],
        "P2_Feature6": [aa_features[1][4]],
        "P1_Degree": [network_features[0][0]],
        "P1_Clustering": [network_features[0][1]],
        "P1_ANN": [network_features[0][2]],
        "P2_Degree": [network_features[1][0]],
        "P2_Clustering": [network_features[1][1]],
        "P2_ANN": [network_features[1][2]],
    }

    df = save_protein_features(data)

    # Predict using trained model
    predict_protein_pairs(
        df=df,
        model_file="./protein_thermo_meso_classifier.pkl"
    )


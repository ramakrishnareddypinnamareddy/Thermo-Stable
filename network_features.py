#!/usr/bin/env python3
import sys
import os
import networkx as nx
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

def calculate_distance(atom1, atom2):
    return atom1 - atom2

def create_protein_network(pdb_file, distance_threshold=6.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    G = nx.Graph()

    for model in structure:
        for chain in model:
            residues = [r for r in chain if is_aa(r) and 'CA' in r]
            for i, res1 in enumerate(residues):
                for res2 in residues[i+1:]:
                    dist = calculate_distance(res1['CA'], res2['CA'])
                    if dist <= distance_threshold:
                        G.add_edge(res1, res2)
    return G

def compute_network_features(pdb_file):
    G = create_protein_network(pdb_file)
    avg_degree = sum(dict(G.degree()).values()) / len(G) if len(G) > 0 else 0
    clustering = nx.average_clustering(G) if len(G) > 0 else 0
    ann_degree = sum(nx.average_neighbor_degree(G).values()) / len(G) if len(G) > 0 else 0
    return avg_degree, clustering, ann_degree

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 network_features.py <protein1.pdb> <protein2.pdb>")
        sys.exit(1)

    pdb1, pdb2 = sys.argv[1], sys.argv[2]

    p1_deg, p1_clust, p1_ann = compute_network_features(pdb1)
    p2_deg, p2_clust, p2_ann = compute_network_features(pdb2)

    # Print in a parseable format for run_classifier
    print(f"{p1_deg}\t{p1_clust}\t{p1_ann}\t{p2_deg}\t{p2_clust}\t{p2_ann}")


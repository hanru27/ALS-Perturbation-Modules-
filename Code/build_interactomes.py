# code for assembling perturbed interactome and random networks as conteols 

import pandas as pd
import random
import os


# Directories 
INTERACTOME_FILE = "Interactome.txt"
PERTURBATION_LIST = "perturbation_list.txt"
OUTPUT_DIR = ""          
RANDOM_NETWORKS_DIR = ""
RANDOM_SEED= 42          

# functions for interactome building and dictionary creation 

def load_interactome(path: str) -> pd.DataFrame:
    df = pd.read_csv(INTERACTOME_FILE, sep="\t", header=None, names=["GeneID_A", "GeneID_B"])
    df["GeneID_A"] = df["GeneID_A"].astype(int)
    df["GeneID_B"] = df["GeneID_B"].astype(int)
    interactome = interactome.drop(columns=2)
    print(f"Loaded interactome: {len(df):,} edges")
    return df


def load_perturbation_list(PERTURBATION_LIST: str) -> set:
    with open(PERTURBATION_LIST) as fh:
        ids = {int(line.strip()) for line in fh if line.strip()}
    print(f"Loaded {len(ids):,} perturbed gene IDs")
    return ids


def build_gene_dictionary(df: pd.DataFrame) -> set:
    all_genes = set(df["GeneID_A"]).union(set(df["GeneID_B"]))
    print(f"Gene universe: {len(all_genes):,} unique gene IDs")
    return all_genes


def remove_genes(df: pd.DataFrame, genes_to_remove: set) -> pd.DataFrame:
    mask = df["GeneID_A"].isin(genes_to_remove) | df["GeneID_B"].isin(genes_to_remove)
    filtered = df[~mask].copy()
    print(f"  Removed {mask.sum():,} edges → {len(filtered):,} edges remaining")
    return filtered

def save_network(df: pd.DataFrame, filename: str) -> None:
    path = os.path.join(OUTPUT_DIR, filename)
    df.to_csv(path, sep="\t", index=False, header=False)
    print(f"  Saved → {path}")


def add_edge_weights(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["weight"] = 1
    return out



###-------------
#Building the networks  

def main():
    rng = random.Random(RANDOM_SEED)

    # Loading the data 
    interactome = load_interactome(INTERACTOME_FILE)
    perturbed_genes = load_perturbation_list(PERTURBATION_FILE)
    gene_dict = build_gene_dictionary(interactome)

   
    regulatory_network = add_edge_weights(interactome)
    save_network(regulatory_network, "Interactome_weighted.txt")

    # Perturbed network = remove all edges touching a perturbed gene (with assurance that all perturbed genes are removed within the network)
    perturbed_df = remove_genes(interactome, perturbed_genes)
    perturbed_df = add_edge_weights(perturbed_df)
    save_network(perturbed_df, "perturbed_interactome.txt")

# checking whether any perturbed genes are still remaining 
    remaining_a = set(perturbed_df["GeneID_A"]) & perturbed_genes
    remaining_b = set(perturbed_df["GeneID_B"]) & perturbed_genes

    if remaining_a or remaining_b:
    print(f"QC FAILED: {len(remaining_a | remaining_b)} perturbed gene(s) still present in perturbed network")
    else:
    print("QC PASSED: no perturbed genes present in perturbed network")




    # Random control networks (generating 5, with assurance that genes from the perturbed list are confidently present within random networks)
   non_perturbed_pool = list(gene_universe - perturbed_genes)
    n_to_remove = len(perturbed_genes)

    for i in range(1, RANDOM_NETWORKS + 1):
        random_genes = set(rng.sample(non_perturbed_pool, n_to_remove))
         mask_r = interactome["GeneID_A"].isin(random_genes) | interactome["GeneID_B"].isin(random_genes)
        random_df = interactome[~mask_r].copy()
        print(f"Random network {i}: removed {mask_r.sum():,} edges → {len(random_df):,} remaining")
        random_df.to_csv(f"random_{i}_interactome.txt", sep="\t", index=False, header=False)
        print(f"Saved → random_{i}_interactome.txt")

    # QC check to verify all perturbed genes are retained in random network
    genes_in_random = set(random_df["GeneID_A"]) | set(random_df["GeneID_B"])
    missing = perturbed_genes - genes_in_random

    if missing:
        print(f"  QC FAILED: {len(missing)} perturbed gene(s) missing from random network {i}")
    else:
        print(f"  QC PASSED: all perturbed genes retained in random network {i}")

print(" All networks saved.")


if __name__ == "__main__":
    main()

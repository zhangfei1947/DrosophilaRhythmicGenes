import os
import csv
import pickle
import lmdb
import pandas as pd

SAMPLES = ["T25LD", "T18DD", "T25DD", "T29DD"]
DATA_DIR = "data"
LMDB_DIR = "lmdbs"

def parse_timepoint_columns(header):
    timepoints = []
    for col in header[1:]:
        if "." in col:
            tp, rep = col.split(".")
            timepoints.append(int(tp))
    return sorted(set(timepoints))

def build_gene_data(csv_path):
    df = pd.read_csv(csv_path, sep="\t", index_col=0)
    timepoints = parse_timepoint_columns(df.columns.tolist())
    gene_dict = {}

    for gene_id, row in df.iterrows():
        expr = {}
        for tp in timepoints:
            r1 = row.get(f"{tp}.r1")
            r2 = row.get(f"{tp}.r2")
            mean_val = (float(r1) + float(r2)) / 2 if pd.notna(r1) and pd.notna(r2) else None
            expr[tp] = {
                "r1": float(r1) if pd.notna(r1) else 0,
                "r2": float(r2) if pd.notna(r2) else 0,
                "mean": mean_val
            }
        gene_dict[gene_id] = expr
    return gene_dict

def write_to_lmdb(gene_data, lmdb_path):
    env = lmdb.open(lmdb_path, map_size=50 * 1024 * 1024)
    with env.begin(write=True) as txn:
        for gene_id, expr in gene_data.items():
            txn.put(gene_id.encode(), pickle.dumps(expr))

if __name__ == "__main__":
    for sample in SAMPLES:
        csv_file = f"{sample}.d4d5.fpkm.csv"
        csv_path = os.path.join(DATA_DIR, csv_file)
        lmdb_path = os.path.join(LMDB_DIR, sample + ".lmdb")
        print(f"[INFO] Processing {csv_file}...")
        gene_data = build_gene_data(csv_path)
        write_to_lmdb(gene_data, lmdb_path)

#################################Making a blast index#########################

import csv
from subprocess import run
import pandas as pd
import os

def make_unique_index(index):
    seen = {}
    unique_index = []
    for idx in index:
        count = seen.get(idx, 0)
        if count:
            unique_idx = f"{idx}-{count + 1}"
            while unique_idx in seen:
                count += 1
                unique_idx = f"{idx}-{count + 1}"
        else:
            unique_idx = idx
        unique_index.append(unique_idx)
        seen[idx] = count + 1
    return unique_index

def df_to_fasta(tab, fasta_filename,id_cols=['Enhancer_ID'],seq_col='Enhancer_sequence'):
    """Takes a pandas df and writes sequences to a FASTA file (give whole path)."""
    tab=tab.loc[~tab[seq_col].isna(),:]
    tab.index=tab[id_cols].astype(str).agg('_'.join, axis=1)
    tab.index = make_unique_index(tab.index)
    with  open(fasta_filename, 'w') as fasta_file:
        for row in tab.index:
            # Adjust indices as per your CSV structure
            fasta_file.write(f'>{row}\n{tab.loc[row,seq_col]}\n')

def create_blast_db(fasta_filename, db_name,env_bin_path=''):
    """Creates a BLAST database from the given FASTA file."""
    cmd = [os.path.join(env_bin_path,'makeblastdb'), '-in', fasta_filename, '-dbtype', 'nucl', '-out', db_name]
    run(cmd)

def query_blast_db(test_seq, db_name, output_file,other_args=[],env_bin_path=''):
    """Queries the BLAST database with the given test sequence. Outputs in commented tabular format (outfmt7)"""
    cmd = [os.path.join(env_bin_path,'blastn'), '-query', test_seq, '-db', db_name, '-out', output_file,'-outfmt', '7']+other_args
    run(cmd)

def read_blast_output7(output):
    """Reads BLAST outfmt 7 to pandas dataframe."""
    return(pd.read_csv(output,comment='#',sep='\t',header=None,
               names=["qseqid", "sseqid", "pident", "length", "mismatch",
                      "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]))

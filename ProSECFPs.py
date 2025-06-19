logo = """
           F H V F V R D V L Q H V D S M Q K D Y P G L P V F  
          _____           _____ ______ _____ ______ _____     
         |  __ \         / ____|  ____/ ____|  ____|  __ \    
         | |__) | __ ___| (___ | |__ | |    | |__  | |__) |__ 
         |  ___/ '__/ _ \\___ \|  __|| |    |  __| |  ___/ __|
         | |   | | | (_) |___) | |___| |____| |    | |   \__ \\
         |_|   |_|  \___/_____/|______\_____|_|    |_|   |___/
                                                                      
          Protein Sequence Extended-Connectivity Fingerprints                                            
        
           0 1 0 0 0 1 1 0 1 0 1 0 1 1 0 0 0 0 1 1 0 0 0 0 0

University of Pisa - Department of Pharmacy - MMVSL https://www.mmvsl.it/wp/
Telethon Institute of Genetics and Medicine - Scuola Superiore Meridionale
\n\n"""

print(logo)

import argparse
parser = argparse.ArgumentParser()

parser.add_argument(
    "-nj", "--n_jobs", 
    type=int, 
    default=1,
    help="Number of parallel jobs to run. Default is 1."
)

parser.add_argument(
    "-in", "--input_dataset", 
    type=str, 
    default="input.csv",
    help="Path to the input CSV file containing preprocessed protein sequences. Default is 'input.csv'."
)

parser.add_argument(
    "-out", "--output_dataset", 
    type=str, 
    default="output.csv",
    help="Path to the output CSV file where the computed descriptors will be saved. Default is 'output.csv'. The output will contain two columns: 'Sequences' (the input protein sequences) and 'C-ProSECFPs' (the corresponding descriptor vectors)."
)

args = parser.parse_args()
nj, input_dataset, output_dataset = args.n_jobs, args.input_dataset, args.output_dataset

import pandas as pd
import numpy as np
import joblib
import hashlib
import itertools
from pandarallel import pandarallel

# Configurable parameters
radius = 12
bits = 1024
sequence_col = 'Sequences'
fingerprint_type = 'Count-MorganFingerprint'  
desc_dataset_path = 'descriptors.dump'

# Initialization of pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=nj)

# Loading dataset
dataset = pd.read_csv(input_dataset)

# Utility functions
def create_fragments(sequence, n):
    return [sequence[max(0, i - n):min(len(sequence), i + n + 1)] for i in range(len(sequence))]

def add_fragments_column(df, seq_col, n):
    df[f'{seq_col}_{n}aa'] = df[seq_col].apply(lambda seq: create_fragments(seq, n))
    return df

def calculate_unaa(sequence):
    return list(sequence)

def hash_function(identifier, n_bits):
    vector_string = ','.join(map(str, identifier))
    sha256 = hashlib.sha256(vector_string.encode('utf-8')).hexdigest()
    return int(sha256, 16) % n_bits

def map_sequence(seq, aa_dict):
    return [aa_dict.get(aa, 'unknown') for aa in seq]

def calculate_hashes(seq, n):
    hashed = []
    for i in range(len(seq)):
        arr = [n, seq[i]]
        arr.extend(itertools.chain.from_iterable([[1, seq[i - j]] for j in range(1, n) if i - j >= 0]))
        arr.extend(itertools.chain.from_iterable([[1, seq[i + j]] for j in range(1, n) if i + j < len(seq)]))
        hashed.append(hash_function(arr, bits))
    return hashed

def trim_list(lst, trim_count):
    return lst[trim_count:-trim_count] if isinstance(lst, list) and len(lst) >= 2 * trim_count else []

def safe_convert_to_int(arr):
    return np.array([int(x) if str(x).isdigit() else 0 for x in arr])

def calculate_count_fingerprint(arr):
    binary_array = np.zeros(bits, dtype=np.int32)
    for num in arr:
        binary_array[int(num) - 1] += 1
    return binary_array

def calculate_unique_fingerprint(arr):
    binary_array = np.zeros(bits, dtype=np.int32)
    for num in set(arr):
        binary_array[int(num) - 1] = 1
    return binary_array

# Initial Preprocessing
dati = pd.DataFrame()
for j in range(1, radius + 1):
    dati = add_fragments_column(dataset, sequence_col, j)

dati[f'{sequence_col}_0aa'] = dataset[sequence_col].parallel_apply(calculate_unaa)

# Loading Normalized Descriptors
df_norm = joblib.load(desc_dataset_path)
hashed_vectors = df_norm.apply(lambda row: hash_function(row.values, bits), axis=1)
df_norm['hash_sha256'] = hashed_vectors
aa_dict = hashed_vectors.to_dict()

# Mapping and Hashing
dati['0_aa_hash'] = dati[f'{sequence_col}_0aa'].parallel_apply(lambda seq: map_sequence(seq, aa_dict))
for i in range(1, radius + 1):
    prev_col = f'{i - 1}_aa_hash'
    dati = dati if prev_col in dati else dati.assign(**{prev_col: dati['0_aa_hash']})
    dati[f'{i}_aa_hash'] = dati[prev_col].parallel_apply(calculate_hashes, args=(i,))

hash_cols = ['0_aa_hash'] + [f'{i}_aa_hash' for i in range(1, radius + 1)]
for idx, col in enumerate(hash_cols[2:], start=1):
    dati[col] = dati[col].apply(lambda lst: trim_list(lst, idx))

dati['Concatenated'] = dati[hash_cols].parallel_apply(lambda row: np.concatenate(row.values), axis=1)
dati['Concatenated'] = dati['Concatenated'].apply(safe_convert_to_int)

# Calculation of Fingerprints
if fingerprint_type == 'Count-MorganFingerprint':
    dati['C-ProSECFPs'] = dati['Concatenated'].parallel_apply(calculate_count_fingerprint)
    output_df = dati[[sequence_col, 'C-ProSECFPs']]
else:
    dati['B-ProSECFPs'] = dati['Concatenated'].parallel_apply(calculate_unique_fingerprint)
    output_df = dati[[sequence_col, 'B-ProSECFPs']]

output_df.to_csv(output_dataset)
print(output_df)
print(f"Output saved as '{output_dataset}'")




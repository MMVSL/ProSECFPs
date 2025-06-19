# ProSECFPs
ProSECFPs is a tool for generating representations of protein sequences. This repository provides a Python script for computing ProSECFP descriptors using a predefined environment with all required dependencies.

# Installation
To set up the required environment, use the provided YAML file:

conda env create -f prosecfps_env.yml
Then, activate the environment:

conda activate prosecfps_env

# Usage
Run the prediction script with Python:

python ProSECFPs.py -in input.csv -out output.csv -nj 1

# Input
The input file must be a CSV containing a column named Sequences, which represents the protein sequences.
A sample CSV file (input.csv) is included in the repository for testing purposes, along with the amino acid descriptor dataset (descriptors.dump) required for computing the ProSECFP representations.

# Output
The script generates a CSV file containing the ProSECFP descriptors for each protein sequence. The descriptors are computed using the C-PSECFP variant, with an iteration radius of 12 and a representation vector length of 1024.

# Dependencies
All necessary dependencies are included in prosecfps_env.yml.

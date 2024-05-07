#! /bin/bash

# This script takes a list of peptides and the associated intensities as input and runs the whole Peptonizer pipeline
# which finally produces a rating of the taxa and their associated scores.

# On first use of the script, run the following to install all dependencies into a new Conda environment
# conda env create -f env.yml
# Then do, conda activate peptonizer

# Step 1: get taxonomy from pout
python unipept_get_taxonomy_from_pout.py --unipept-response-file "../data/unipept_response.json" --pep-out "../data/peptides.out" --taxonomy-query "2" --fdr "0.01" --pout-file "../sample_input/rescored.psms.tsv" --log-file "../data/unipept.log"

# Step 2: compute weighted taxa
python weight_taxa.py --unipept-response-file "../data/unipept_response.json" --number-of-taxa "50" --out "../data/weights_data_frame.csv" --unipept-peptides "../data/unipept_peptides.json" --taxa-rank "species"

# Step 3: run the believe

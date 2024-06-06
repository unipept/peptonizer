#! /bin/bash

# This script takes a list of peptides and the associated intensities as input and runs the whole Peptonizer pipeline
# which finally produces a rating of the taxa and their associated scores.

# On first use of the script, run the following to install all dependencies into a new Conda environment
# Make sure to start from a Docker image where mamba is already installed!
# 1) conda env create -f env.yml (make sure to run this command from the workflow directory)
# 2) mamba install -c conda-forge -c bioconda -n peptonizer snakemake
# 3) conda activate peptonizer
# 4) snakemake --use-conda --cores 1

# The input file that this script requires is a pout-file with rescored peptides
# By default, it is called "rescored.psms.tsv"

# Step 1: get taxonomy from pout
python unipept_get_taxonomy_from_pout.py --unipept-response-file "../../data/unipept_response.json" --unipept-peptide-counts "../../data/unipept_peptide_counts.json" --taxonomy-query "2" --fdr "0.01" --pout-file "../../sample_input/rescored.psms.tsv" --log-file "../../data/unipept.log"

# Step 2: compute weighted taxa
python weight_taxa.py --unipept-response-file "../../data/unipept_response.json" --number-of-taxa "50" --out "../../data/weights_data_frame.csv" --taxa-weight-file "../../data/taxa_weights.csv"  --unipept-peptide-counts "../../data/unipept_peptide_counts.json" --taxon-rank "species"

# Step 3: create a factor graph using the weights that have been calculated in the previous step.
python3 create_pepgm_graph.py --graph-data-frame "../../data/weights_data_frame.csv" --out "../../data/pepgm_graph.graphml"

# Step 4: actually run the PepGM algorithm itself
python3 pepgm.py --graphml-path "../../data/pepgm_graph.graphml" --prior "0.1" --alpha "0.7" --beta "0.2" --regularized "True" --out "../../data/pepgm_result.csv"

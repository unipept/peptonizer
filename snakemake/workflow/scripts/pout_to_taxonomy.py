import argparse
import gzip
import json

from peptonizer.peptonizer import parse_ms2rescore_output, fetch_unipept_taxon_information


parser = argparse.ArgumentParser()

parser.add_argument(
    "--taxonomy-query",
    required=True,
    help="Taxa that should be used to query in Unipept. If querying all taxa, put [1].",
)
parser.add_argument(
    "--fdr",
    type=float,
    required=True,
    help="Min peptide score for the peptide to be included in the search.",
)
parser.add_argument(
    "--pout-file",
    type=str,
    required=True,
    help="Input: path to percolator (ms2rescore) '.pout' file.",
)
parser.add_argument(
    "--unipept-response-file",
    type=str,
    required=True,
    help="Path to output file that contains all queried peptide counts (which should be used in the next step)."
)
parser.add_argument(
    "--log-file",
    type=str,
    required=True,
    help="Output: path to logfile where failed Unipept query attempts are stored.",
)
parser.add_argument(
    "--taxon-rank",
    type=str,
    required=False,
    default="species",
    help="Taxonomic rank at which you want the Peptonizer2000 results to be resolved.",
)

args = parser.parse_args()

with gzip.open(args.pout_file, 'rt', encoding='utf-8') as file:
    file_contents = file.read()

# Parse the input MS2Rescore file
pep_score, pep_psm_counts = parse_ms2rescore_output(file_contents, args.fdr)

unipept_response = fetch_unipept_taxon_information(
    list(pep_score.keys()),
    args.taxonomy_query,
    args.taxon_rank,
    args.log_file
)

with open(args.unipept_response_file, "w") as f:
    f.write(json.dumps(unipept_response))

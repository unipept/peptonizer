import argparse

from peptonizer import parse_pout, parse_ms2rescore_output, fetch_unipept_taxon_information


parser = argparse.ArgumentParser()

parser.add_argument(
    "--unipept-response-file",
    type=str,
    required=True,
    help="Output: path to Unipept response .json file",
)
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
    nargs="+",
    required=True,
    help="Input: paths to percolator (ms2rescore) '.pout' files.",
)
parser.add_argument(
    "--unipept-peptide-counts",
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

args = parser.parse_args()

pep_score_psm = parse_ms2rescore_output(args.pout_file, args.fdr, "")
fetch_unipept_taxon_information(
    pep_score_psm,
    args.unipept_peptide_counts,
    args.unipept_response_file,
    args.taxonomy_query,
    args.log_file
)

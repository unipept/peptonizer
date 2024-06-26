import argparse

from peptonizer import perform_taxa_weighing


parser = argparse.ArgumentParser()

parser.add_argument(
    "--unipept-response-file",
    type=str,
    required=True,
    help="Input: path to a Unipept response '.json' file that's been produced earlier in the pipeline.",
)
parser.add_argument(
    "--number-of-taxa",
    type=int,
    required=True,
    help="Number of taxa to include in the final Peptonizer2000 output.",
)
parser.add_argument("--out", type=str, required=True, help="path to csv out file")
parser.add_argument(
    "--taxa-weight-file",
    type=str,
    required=False,
    help="Output: path to a CSV-file that will contain all computed taxa weights.",
)
parser.add_argument(
    "--unipept-peptide-counts",
    type=str,
    required=True,
)
parser.add_argument(
    "--taxon-rank",
    type=str,
    required=False,
    default="species",
    help="Taxonomic rank at which you want the Peptonizer2000 results to be resolved.",
)

args = parser.parse_args()

df, weights = perform_taxa_weighing(
    args.unipept_response_file,
    args.unipept_peptide_counts,
    args.number_of_taxa,
    taxa_rank=args.taxon_rank,
)

print("Started dumping produced results to CSV-files...")
df.to_csv(args.out)
weights.to_csv(args.taxa_weight_file)

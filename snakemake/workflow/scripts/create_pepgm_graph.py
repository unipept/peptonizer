import argparse
import pandas as pd

from peptonizer.peptonizer import generate_pepgm_graph


parser = argparse.ArgumentParser(
    description="Run the PepGM algorithm from command line"
)

parser.add_argument(
    "--sequence-scores-dataframe-file",
    type=str,
    required=True,
    help="Dataframe file containing the taxa weights that have been computed before.",
)
parser.add_argument(
    "--out",
    type=str,
    required=True,
    help="Path to output file where GraphML will be saved.",
)
parser.add_argument(
    "--use-communities",
    type=bool,
    required=False,
    default=True
)

args = parser.parse_args()

ct_factor_graph = generate_pepgm_graph(pd.read_csv(
    args.sequence_scores_dataframe_file),
    use_communities=args.use_communities
)
ct_factor_graph.save_to_graph_ml(args.out)

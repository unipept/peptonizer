import argparse

from peptonizer import generate_pepgm_graph


parser = argparse.ArgumentParser(
    description="Run the PepGM algorithm from command line"
)

parser.add_argument(
    "--graph-data-frame",
    type=str,
    required=True,
    help="Path to where you want to save the GraphML file of the factor graph.",
)
parser.add_argument(
    "--out",
    type=str,
    required=True,
    help="Path to output file where GraphML will be saved.",
)

args = parser.parse_args()

generate_pepgm_graph(args.graph_data_frame, args.out)

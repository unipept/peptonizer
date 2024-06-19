import argparse
from factor_graph_generation import *

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

taxon_graph = TaxonGraph()
taxon_graph.create_from_unipept_response_csv(args.graph_data_frame)
factor_graph = FactorGraph()
factor_graph.construct_from_existing_graph(taxon_graph)
ct_factor_graph = generate_ct_factor_graphs(factor_graph)
ct_factor_graph.save_to_graph_ml(args.out)

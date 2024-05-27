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
taxon_graph.CreateFromUnipeptResponseCSV(args.graph_data_frame)
factor_graph = FactorGraph()
factor_graph.ConstructFromTaxonGraph(taxon_graph)
ct_factor_graph = GenerateCTFactorGraphs(factor_graph)
ct_factor_graph.SaveToGraphML(args.out)

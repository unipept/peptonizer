import argparse
import networkx as nx

from zero_lookahead_belief_propagation import *
from factor_graph_generation import *


parser = argparse.ArgumentParser(
    description="Run the PepGM algorithm from command line"
)

parser.add_argument(
    "--graphml-path",
    type=str,
    required=True,
    help="Path to where the GraphML file of the factor graph is stored.",
)
parser.add_argument(
    "--max-iter",
    nargs="?",
    type=int,
    default=10000,
    help="Max. number of iterations the belief propagation algo will go through.",
)
parser.add_argument(
    "--tol",
    nargs="?",
    type=float,
    default=0.006,
    help="Residual error allowed for the BP algorithm.",
)
parser.add_argument(
    "--out",
    type=str,
    required=True,
    help="Path to the file you want to save your results as.",
)
parser.add_argument(
    "--alpha",
    type=float,
    required=True,
    help="Detection probability of a peptide for the noisy-OR model.",
)
parser.add_argument(
    "--beta", type=float, required=True, help="Probability of wrong detection."
)
parser.add_argument(
    "--prior", type=float, required=True, help="Prior assigned to all taxa."
)
parser.add_argument(
    "--regularized",
    type=bool,
    default=False,
    help="If True, the regularized version of the noisy-OR model is used.",
)

args = parser.parse_args()

ct_factor_graph = CTFactorGraph(args.graphml_path)
ct_factor_graph.FillInFactors(args.alpha, args.beta, args.regularized)
ct_factor_graph.FillInPriors(args.prior)
ct_factor_graph.AddCTNodes()


ct_factor_graphs = [
    SeparateSubgraphs(ct_factor_graph, filter_nodes)
    for filter_nodes in nx.connected_components(ct_factor_graph)
]

results_dict, node_types = calibrate_all_subgraphs(
    ct_factor_graphs, args.max_iter, args.tol
)
save = save_results_to_csv(results_dict, node_types, args.out)

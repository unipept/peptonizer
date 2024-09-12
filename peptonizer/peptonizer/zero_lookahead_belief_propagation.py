# implementation of belief propagation on a peptide-protein graph
# __________________________________________________________________________________________
import time
from sys import getsizeof

from .convolution_tree import *
from .factor_graph_generation import *
from .pqdict import pqdict
from typing import Dict, Any, List, Tuple, Iterator


class Messages:
    """
    Class holding all messages and beliefs.
    Functions execute loopy residual belief propagation with zero look ahead
    """

    # class that holds the messages of iteration t and iteration t+1 as dictionaries
    def __init__(self, ct_graph_in: CTFactorGraph):
        # TODO check if I truly need all three of msg new, msglog and msg.
        self.max_val: Optional[Tuple[int, int]] = None
        self.priorities: pqdict = pqdict({}, reverse=True)
        self.category: str = ct_graph_in.category

        # Maps a node identifier (as specified by the CTGraph) onto a unique integer.
        self.node_descriptions: Dict[str, int] = {}
        # Reverse mapping of the dict above
        self.node_id_to_description: List[str] = []
        # Maps a node ID onto its category
        self.categories: List[str] = []
        # Maps a node ID onto a list of its neighbouring node IDs
        self.neighbours: List[List[int]] = []
        # Keeps track of the number of parents of a node
        self.number_of_parents: List[int] = []
        # Keep track of all edges (with node IDs)
        self.edges: List[Tuple[int, int]] = []

        # Keeps track of residuals for edges in the graph
        self.full_residuals: Dict[Tuple[int, int], float] = {}
        self.total_residuals: Dict[Tuple[Tuple[int, int], Tuple[int, int]], float] = {}

        # Maps a node ID onto its initial belief value
        self.initial_beliefs: List[npt.NDArray[np.float64]] = []
        # Maps a node ID onto its current belief value
        self.current_beliefs: List[npt.NDArray[np.float64]] = []

        self.msg: Dict[Tuple[int, int], npt.NDArray[np.float64]] = {}
        self.msg_new: Dict[Tuple[int, int], npt.NDArray[np.float64]] = {}
        self.msg_log: Dict[Tuple[int, int], npt.NDArray[np.float64]] = {}

        nodes: Iterator[Tuple[Any, Dict[str, Any]]] = ct_graph_in.nodes(data=True)

        for (node_id, node) in enumerate(nodes):
            # Each node should occur exactly once
            assert node[0] not in self.node_descriptions

            self.node_descriptions[node[0]] = node_id
            self.node_id_to_description.append(node[0])
            self.categories.append(node[1]["category"])

            # Only convolution trees have a number of parents property assigned to them.
            if node[1]["category"] == "convolution_tree":
                self.number_of_parents.append(node[1]["NumberOfParents"])
            else:
                self.number_of_parents.append(0)

            if node[1]["category"] == "factor":
                self.initial_beliefs.append(node[1]["InitialBelief"].factor.array)
            elif node[1]["category"] == "peptide" or node[1]["category"] == ct_graph_in.category:
                self.initial_beliefs.append(np.array([node[1]["InitialBelief_0"], node[1]["InitialBelief_1"]]))
            else:
                # this entry will never be used as the convolution trees do not hold beliefs
                self.initial_beliefs.append(np.ones(4))

        # We start out with current beliefs that are identical to the initial beliefs
        self.current_beliefs = self.initial_beliefs.copy()

        # Now that all nodes have been processed, we need to replace all the neighbours of a node by their node IDs
        for (node_id, node) in enumerate(ct_graph_in.nodes(data=True)):
            self.neighbours.append([self.node_descriptions[n] for n in ct_graph_in.neighbors(node[0])])

        # Now, also replace the edge descriptions by the corresponding node IDs
        for start_node, end_node, data in ct_graph_in.edges(data=True):
            start_node_id = self.node_descriptions[start_node]
            end_node_id = self.node_descriptions[end_node]

            self.edges.append((start_node_id, end_node_id))

            if "MessageLength" in data:
                self.msg[(start_node_id, end_node_id)] = np.ones(data["MessageLength"])
                self.msg_new[(start_node_id, end_node_id)] = np.ones(data["MessageLength"])
            else:
                self.msg[(start_node_id, end_node_id)] = np.array([0.5, 0.5])
                self.msg_new[(start_node_id, end_node_id)] = np.array([0, 0])

            self.msg[(end_node_id, start_node_id)] = self.msg[(start_node_id, end_node_id)]
            self.msg_new[(end_node_id, start_node_id)] = self.msg_new[(start_node_id, end_node_id)]

        self.msg_log = self.msg_new.copy()

    # variables (peptides,proteins,taxa)
    def get_incoming_message_variable(self, node: int, node_in: int) -> npt.NDArray[np.float64]:
        return self.msg[node, node_in]

    def compute_out_message_variable(self, node_out: int, node_in: int) -> npt.NDArray[np.float64]:
        incoming_messages: List[npt.NDArray[np.float64]] = []
        node_belief: npt.NDArray[np.float64] = self.current_beliefs[node_out]

        for node_out_neighbor in self.neighbours[node_out]:
            if node_out_neighbor != node_in:
                incoming_messages.append(
                    self.get_incoming_message_variable(node_out_neighbor, node_out)
                )

        if not incoming_messages:
            # TODO Tanja: what exactly do you want to compare here?
            return node_belief if any(node_belief == self.initial_beliefs[node_out]) else self.msg[node_out, node_in]

        # need for logs to prevent underflow in very large multiplications
        incoming_messages_array = np.asarray(np.log(incoming_messages)).reshape(
            len(incoming_messages), 2
        )

        out_message_log = array_utils.log_normalize(
            np.asarray(
                [
                    np.sum([np.log(node_belief[0]), np.sum(incoming_messages_array[:, 0])]),
                    np.sum([np.log(node_belief[1]), np.sum(incoming_messages_array[:, 1])]),
                ]
            )
        )

        if not np.all(out_message_log):
            out_message_log[out_message_log == 0] = 1e-30

        return out_message_log

    # factors (Conditional probability tables), handles different dimension of output/input variables
    def get_incoming_message_factor(self, node: int, node_in: int) -> npt.NDArray[np.float64]:
        return self.msg[node, node_in]

    def compute_out_message_factor(self, node_out: int, node_in: int) -> npt.NDArray[np.float64]:
        incoming_messages: List[npt.NDArray[np.float64]] = []
        node_belief = self.current_beliefs[node_out]

        for node_out_neighbour in self.neighbours[node_out]:
            if node_out_neighbour != node_in:
                if [self.get_incoming_message_factor(node_out_neighbour, node_out)]:
                    # only the messages that have changed get multiplied into the current belief again
                    incoming_messages.append(
                        self.get_incoming_message_factor(node_out_neighbour, node_out)
                    )

        if self.categories[node_in] == "convolution_tree":
            # handles empty & messages with only one value
            incoming_messages.append(np.asarray([1.0, 1.0]))
            in_messages_array: npt.NDArray[np.float64] = np.asarray(incoming_messages).reshape(
                len(incoming_messages), 2
            )  # np.asarray(np.log(IncomingMessages)).reshape(len(IncomingMessages),2)
            out_messages = array_utils.normalize(
                np.multiply(
                    node_belief,
                    [np.prod(in_messages_array[:, 0]), np.prod(in_messages_array[:, 1])],
                )
            )  # lognormalize(np.add(np.log(NodeBelief),[np.sum(IncomingMessages[:,0]),np.sum(IncomingMessages[:,1])]))#

            return np.add(out_messages[:, 0], out_messages[:, 1])
        else:
            if np.asarray(incoming_messages[0]).shape[0] > 2:
                incoming_messages_log = np.log(
                    np.asarray(incoming_messages).reshape(
                        incoming_messages[0].shape[0], 1
                    )
                )
                # catch warning for log(0)

                log_belief = np.log(node_belief)
                out_messages_log = array_utils.log_normalize(np.add(log_belief, incoming_messages_log))
                if not np.all(out_messages_log):
                    out_messages_log[out_messages_log == 0] = 1e-30

                return np.asarray([np.sum(out_messages_log[0, :]), np.sum(out_messages_log[1, :])])
            else:
                incoming_messages.append(np.asarray([1.0, 1.0]))
                incoming_messages_array = np.asarray(incoming_messages).reshape(
                    len(incoming_messages), 2
                )
                out_messages = array_utils.normalize(
                    np.multiply(
                        node_belief,
                        [
                            np.prod(incoming_messages_array[:, 0]),
                            np.prod(incoming_messages_array[:, 1]),
                        ],
                    )
                )
                return np.asarray([np.sum(out_messages[0, :]), np.sum(out_messages[1, :])])

    def compute_out_messages_ct_tree(self, node: int):
        prot_prob_list: List[npt.NDArray[np.float64]] = []
        old_prot_prob_list: List[npt.NDArray[np.float64]] = []
        shared_likelihoods: npt.NDArray[np.float64] = np.ones(self.number_of_parents[node] + 1)
        old_shared_likelihoods: npt.NDArray[np.float64] = np.empty(1)
        peptides: List[int] = []
        prot_list: List[int] = []

        for node_in in self.neighbours[node]:
            if self.categories[node_in] != "factor":
                prot_prob_list.append(self.msg[node_in, node])
                old_prot_prob_list.append(self.msg_log[node_in, node])
                prot_list.append(node_in)
            else:
                peptides.append(node_in)
                try:
                    shared_likelihoods = np.multiply(
                        array_utils.avoid_underflow(shared_likelihoods), self.msg[node_in, node]
                    )
                except:
                    print(shared_likelihoods, self.msg[node_in, node])
                old_shared_likelihoods = np.multiply(
                    array_utils.avoid_underflow(shared_likelihoods), self.msg_log[node_in, node]
                )

        if (
                not np.array_equal(old_shared_likelihoods, shared_likelihoods) and
                any(
                    prot_prob_list[i][0] != old_prot_prob_list[i][0]
                    for i in range(len(prot_prob_list))
                )
        ):
            # only update when the shared likelihoods or at least one of the protein messages has changed
            ct = ConvolutionTree(shared_likelihoods, prot_prob_list)

            for protein in range(len(prot_list)):
                self.msg_new[node, prot_list[protein]] = ct.message_to_variable(protein)
                if not np.all(self.msg_new[node, prot_list[protein]]):
                    self.msg_new[node, prot_list[protein]][
                        self.msg_new[node, prot_list[protein]] == 0
                        ] = 1e-30

            for pep in peptides:
                self.msg_new[node, pep] = ct.message_to_shared_likelihood()
                if not np.all(self.msg_new[node, pep]):
                    self.msg_new[node, pep][self.msg_new[node, pep] < 1e-30] = 1e-30

        else:
            for protein in range(len(prot_list)):
                self.msg_new[node, prot_list[protein]] = self.msg[node, prot_list[protein]]

            for pep in peptides:
                self.msg_new[node, pep] = self.msg[node, pep]

    def compute_infinity_norm_residual(self, start_node: int, end_node: int) -> float:
        msg1 = self.msg[start_node, end_node]
        msg2 = self.msg_log[start_node, end_node]

        if len(self.msg_log[start_node, end_node]) != len(self.msg[start_node, end_node]):
            msg2 = [1] * len(self.msg[start_node, end_node])

        pos = 0
        for i in msg1:
            if i < 1e-30:
                msg1[pos] = 1e-30
            pos += 1

        return np.max(np.abs(np.log(np.divide(msg1, msg2))))
        # approximate residual with zero look-ahead

    def compute_zero_look_ahead_residual(self, start_node: int, end_node: int) -> float:
        node_in_neighbours = [nodes for nodes in self.neighbours[start_node]]
        node_in_neighbours.remove(end_node)
        return sum(
            [self.full_residuals[(neighbour, start_node)] for neighbour in node_in_neighbours]
        )

    def compute_total_residuals(self, start_node: int, end_node: int, current_residual: float):
        for start_neighbour in self.neighbours[start_node]:
            if start_neighbour != end_node:
                self.total_residuals[
                    ((start_neighbour, start_node), (start_node, end_node))
                ] = 0

        for end_neighbour in self.neighbours[end_node]:
            if end_neighbour != start_node:
                self.total_residuals[((start_node, end_node), (end_node, end_neighbour))] = (
                    self.total_residuals[((start_node, end_node), (end_node, end_neighbour))] + current_residual
                )

    def compute_priority(self, start_node: int, end_node: int):
        self.priorities[start_node, end_node] = 0

        for end_neighbor in self.neighbours[end_node]:
            if end_neighbor != start_node:
                self.priorities[end_node, end_neighbor] = np.sum(
                    [
                        self.total_residuals[(sum_run, end_node), (end_node, end_neighbor)]
                        for sum_run in self.neighbours[end_node]
                        if sum_run != end_neighbor
                    ]
                )

    # computes new message for a given edge (startname, endname) in the direction startname -> endname
    def single_edge_direction_update(self, start_node: int, end_node: int, checked_cts: set[int]):
        if (
                self.categories[start_node] == self.category
                or self.categories[start_node] == "peptide"
        ):
            self.msg_new[start_node, end_node] = self.compute_out_message_variable(
                start_node, end_node
            )

        if self.categories[start_node] == "convolution_tree" and start_node not in checked_cts:
            self.compute_out_messages_ct_tree(start_node)
            checked_cts.add(start_node)

        if self.categories[start_node] == "factor":
            self.msg_new[start_node, end_node] = self.compute_out_message_factor(
                start_node, end_node
            )

    # compute updated messages for all edges
    def compute_update(self, local_loops: bool = False):
        checked_cts: set[int] = set()  # keeps track of which CT has already been active

        if local_loops and self.max_val:
            for end_node in self.neighbours[self.max_val[1]]:
                start_node = self.max_val[1]
                self.single_edge_direction_update(start_node, end_node, checked_cts)

        else:
            for edge in self.edges:
                # update all edges
                start_node: int = edge[0]
                end_node: int = edge[1]

                self.single_edge_direction_update(start_node, end_node, checked_cts)

                start_node: int = edge[1]
                end_node: int = edge[0]
                self.single_edge_direction_update(start_node, end_node, checked_cts)

    def get_priority_message(self, priority_vector: pqdict) -> Tuple[int, int]:
        self.max_val = priority_vector.top()
        return self.max_val

    def zero_look_ahead_loopy_loop(self, max_loops: int, tolerance: float, local: bool = False) -> Dict[str, npt.NDArray[np.float64]]:
        """
        Run the zero-look-ahead belief propagation algorithm.
        :param max_loops: int, maximum number of iterations in case of non-convergence
        :param tolerance: float, tolerance for convergence check
        :param local: Bool, parameter passed to Computed Update function
        """

        max_residual = 100

        # first, do 5 loops where I update all messages
        print(f"Time spent in loop 0/{max_loops}: 0s")
        for k in range(0, 5):
            start_t = time.time()
            self.compute_update()
            self.msg_log.update(self.msg)
            self.msg.update(self.msg_new)
            k += 1
            end_t = time.time()
            print(f"\rTime spent in loop {k}/{max_loops}: {(end_t - start_t):.3f}s")

        # compute all residuals after 5 runs once (= initialize the residual/priorities vectors)
        for edge in self.edges:
            # compute all residuals of the messages in this loop
            for start_node, end_node in [(edge[1], edge[0]), (edge[0], edge[1])]:
                self.full_residuals[
                    (start_node, end_node)
                ] = self.compute_infinity_norm_residual(start_node, end_node)

                # initialize the total residual to 0
                for end2 in self.neighbours[end_node]:
                    self.total_residuals[((start_node, end_node), (end_node, end2))] = 0
                    self.total_residuals[((end2, end_node), (end_node, start_node))] = 0

            if len(self.total_residuals) % 1000 < 5:
                print(f"Total residuals length: {len(self.total_residuals)}, size in bytes: {getsizeof(self.total_residuals)}")

        # set the priority vector once with copy of the previously calculated residuals
        self.priorities = pqdict(self.full_residuals.copy(), reverse=True)

        k = 5

        print(f"\rTime spent in loop 0/{max_loops}: 0s -> residual max 0", end="")
        while k < max_loops and max_residual > tolerance:
            checked_cts: set[int] = set()
            # actual zero-look-ahead-BP part
            start_t = time.time()
            priority_message = self.get_priority_message(self.priorities)
            max_residual = self.priorities[priority_message]
            self.single_edge_direction_update(priority_message[0], priority_message[1], checked_cts)
            priority_residual = self.compute_infinity_norm_residual(
                priority_message[0], priority_message[1]
            )
            self.msg_log.update(self.msg)
            self.msg.update(self.msg_new)
            self.compute_total_residuals(
                priority_message[0], priority_message[1], priority_residual
            )
            self.compute_priority(priority_message[0], priority_message[1])

            end_t = time.time()

            # Only update the time per loop every 5 iterations
            if k % 5 == 0:
                print(
                    f"\rTime spent in loop {k}/{max_loops}: {(end_t - start_t):.3f}s -> residual max {max_residual:.3f}",
                    end="")

            k += 1

        print()

        # marginalize once the model has converged
        for (node_id, node_category) in enumerate(self.categories):
            if (
                    node_category == self.category
                    or node_category == "peptide"
            ):
                incoming_messages: List[npt.NDArray[np.float64]] = []

                for variable_neighbour in self.neighbours[node_id]:
                    incoming_messages.append(
                        self.get_incoming_message_variable(variable_neighbour, node_id)
                    )

                # log to avoid overflow
                incoming_messages_array = np.asarray(np.log(incoming_messages)).reshape(
                    len(incoming_messages), 2
                )
                logged_variable_marginal = array_utils.log_normalize(
                    np.asarray(
                        [
                            np.sum(
                                [
                                    np.log(self.initial_beliefs[node_id][0]),
                                    np.sum(incoming_messages_array[:, 0]),
                                ]
                            ),
                            np.sum(
                                [
                                    np.log(self.initial_beliefs[node_id][1]),
                                    np.sum(incoming_messages_array[:, 1]),
                                ]
                            ),
                        ]
                    )
                )

                self.current_beliefs[node_id] = logged_variable_marginal

        # Translate the node_ids back to their original sequences and return the current beliefs as a dictionary
        output_beliefs: Dict[str, npt.NDArray[np.float64]] = {}
        for (node_id, beliefs) in enumerate(self.current_beliefs):
            output_beliefs[self.node_id_to_description[node_id]] = beliefs
        return output_beliefs


# calibration through message passing of all subgraphs in the List of factor graphs
def calibrate_all_subgraphs(list_of_ct_factor_graphs: List[CTFactorGraph], max_iterations: int, tolerance: float, local: bool = False) -> Tuple[Dict[str, npt.NDArray[np.float64]], Dict[str, str]]:
    """
    Performs bayesian inference through loopy belief propagation, returns dictionary {variable:posterior_probability}
    :param list_of_ct_factor_graphs: list, contains FactorGraph objects on which inference can be performed
    :param max_iterations: int, max number of iterations in case of non-convergence
    :param tolerance: float, error tolerance between messages for convergence criterion
    :param local: Bool, whether loops are calculated locally
    """

    results_dict: Dict[str, npt.NDArray[np.float64]] = {}
    node_category_dict: Dict[str, str] = {}

    for (idx, graph) in enumerate(list_of_ct_factor_graphs):
        print(f"Started calibrating graph {idx + 1} of {len(list_of_ct_factor_graphs)}")
        if graph.number_of_nodes() > 2:
            node_category_dict.update(dict(graph.nodes(data="category")))
            initialized_message_object = Messages(graph)
            current_beliefs = initialized_message_object.zero_look_ahead_loopy_loop(
                max_iterations, tolerance, local
            )
            results_dict.update(current_beliefs)
    return results_dict, node_category_dict


def convert_results_to_csv(results_dict: Dict[str, npt.NDArray[np.float64]], node_dict: Dict[str, str]):
    """
    Save Loopy Belief Propagation results to .csv file
    :param results_dict: dict, {variable:posterior_probability}
    :param node_dict: dict, dictionary of nodes that were in the factor graph and their attributes, to include the node category in the results
    :param name_string: str, csv output path
    """

    full_results_dict = {
        key: [results_dict[key][1], node_dict[key]] for key in results_dict.keys()
    }

    return pd.DataFrame.from_dict(data=full_results_dict, orient="index").to_csv(header=False)


def run_belief_propagation(
        graphml_content: str,
        alpha: float,
        beta: float,
        regularized: bool,
        prior: float,
        max_iter: int = 10000,
        tol: float = 0.006
    ):
    """
    Runs the belief propagation algorithm on a graph that's represented by the string in graphml_content with the
    tuning parameters further specified to this function. This function returns a string that contains the result of
    the belief propagation algorithm, represented as a CSV (and can thus directly be written to a CSV-file, if desired).
    """

    ct_factor_graph = CTFactorGraph(graphml_content)
    ct_factor_graph.fill_in_factors(alpha, beta, regularized)
    ct_factor_graph.fill_in_priors(prior)
    ct_factor_graph.add_ct_nodes()

    ct_factor_graphs = [
        separate_subgraphs(ct_factor_graph, filter_nodes)
        for filter_nodes in nx.connected_components(ct_factor_graph)
    ]

    results_dict, node_types = calibrate_all_subgraphs(
        ct_factor_graphs, max_iter, tol
    )

    return convert_results_to_csv(results_dict, node_types)

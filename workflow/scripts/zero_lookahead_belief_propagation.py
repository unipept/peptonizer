# implementation of belief propagation on a peptide-protein graph
# __________________________________________________________________________________________
import numpy as np
import math
import json
import pandas as pd

from scipy.signal import fftconvolve
from factor_graph_generation import *
from scipy.special import logsumexp
from pqdict import pqdict
from numba import jit

import time


def normalize(array):
    return array / np.sum(array)


# normalization of log probabilities
def log_normalize(array):
    try:
        y = np.exp(array - logsumexp(array))

    except FloatingPointError:
        print(array)
        y = np.exp(array - logsumexp(array))
    return y


def avoid_underflow(array):
    array[array < 1e-30] = 1e-30
    return array


# implementation of the convolution tree according to serang
# not written by me!!
class CTNode:
    def __init__(self, jointAbove):
        # normalize for greater precision
        self.jointAbove = normalize(jointAbove)

        self.leftParent = None
        self.rightParent = None

        self.likelihoodBelow = None

    # passing msges down: adding variables
    @classmethod
    def createCountNode(cls, lhs, rhs):
        # create a cound node with joint prob for two parents above(vs the init if we have no parents)
        jointAbove = fftconvolve(lhs.jointAbove, rhs.jointAbove)
        result = cls(jointAbove)

        result.leftParent = lhs
        result.rightParent = rhs

        return result

    # passing messages up : subtracting variables
    def messageUp(self, answerSize, otherJointVector):
        startingPoint = len(otherJointVector) - 1
        result = fftconvolve(otherJointVector[::-1], self.likelihoodBelow)[
                 startingPoint: startingPoint + answerSize
                 ]

        return normalize(result)

    def messageUpLeft(self):
        return self.messageUp(
            len(self.leftParent.jointAbove[0]), self.rightParent.jointAbove[0]
        )

    def messageUpRight(self):
        return self.messageUp(
            len(self.rightParent.jointAbove[0]), self.leftParent.jointAbove[0]
        )

    # once all messages are received
    def posterior(self):
        return normalize(self.jointAbove * self.likelihoodBelow)

    def MessagesUp(self):
        check1 = self.jointAbove
        check2 = self.likelihoodBelow
        return self.likelihoodBelow


class ConvolutionTree:
    def __init__(self, nToSharedLikelihoods, proteins):
        self.nToSharedLikelihoods = nToSharedLikelihoods
        self.logLength = int(math.ceil(np.log2(float(len(proteins)))))  # length we need
        self.allLayers = []
        self.buildFirstLayer(proteins)
        self.buildRemainingLayers()
        self.propagateBackward()
        self.nProteins = len(proteins)

    def buildFirstLayer(self, proteins):
        # construct first layer (of proteins)
        layer = []
        for prot in proteins:
            protNode = CTNode(prot)
            layer.append(protNode)

        # pad with necessarily absent dummy variables so that the
        # number of variables is a power of 2; this is not the most
        # efficient method for this. because they are absent, they won't influence the
        # total sum, and thus Ds.
        for i in range(0, 2 ** self.logLength - len(proteins)):
            # this protein cannot be present, therefor set propbaility array to (0,1)
            layer.append(CTNode([np.array([1, 0])]))  # TODO change this order

        self.allLayers.append(layer)

    def buildRemainingLayers(self):
        # construct layers of count nodes
        for L in range(self.logLength):
            # print('layers needed: ',int(len(self.allLayers[0])/(2**(L+1))))
            mostRecentLayer = self.allLayers[-1]
            layer = []
            for i in range(int(len(self.allLayers[0]) / (2 ** (L + 1)))):
                leftParent = mostRecentLayer[i * 2]
                rightParent = mostRecentLayer[i * 2 + 1]
                countNode = CTNode.createCountNode(leftParent, rightParent)
                layer.append(countNode)

            # add connection to remaining nodes (when layer above is not a power of 2)
            self.allLayers.append(layer)

        # final node gets (Ds | N) multiplied into its likelihoodBelow
        finalNode = self.allLayers[-1][0]
        # normalize for greater precision
        finalNode.likelihoodBelow = normalize(self.nToSharedLikelihoods)
        self.LastNode = finalNode

    def propagateBackward(self):
        # propagate backward, setting likelihoodBelow.
        # the loop has upper bound at logLength+1
        # because of the layer of proteins
        for L in range(1, self.logLength + 1)[::-1]:
            layer = self.allLayers[L]

            for i in range(len(layer)):
                node = layer[i]
                leftParent = node.leftParent
                rightParent = node.rightParent

                leftParent.likelihoodBelow = node.messageUpLeft()
                rightParent.likelihoodBelow = node.messageUpRight()

        self.proteinLayer = self.allLayers[0]

    def posteriorForVariable(self, protInd):
        return self.proteinLayer[protInd].posterior()

    def MessageToVariable(self, protInd):
        return self.proteinLayer[protInd].MessagesUp()

    def MessageToSharedLikelihood(self):
        return self.LastNode.jointAbove[0][0: (self.nProteins + 1)]


class Messages:
    """
    Class holding all messages and beliefs.
    Functions execute loopy residual belief propagation with zero look ahead
    """

    # class that holds the messages of iteration t and iteration t+1 as dictionaries

    def __init__(self, ct_graph_in):
        self.msg = {}
        self.msg_new = {}
        self.msg_log = {}
        self.graph = ct_graph_in
        self.max_val = None
        self.full_residual = {}
        self.initial_beliefs = {}
        self.current_beliefs = {}
        self.queue = {}
        self.priorities = pqdict({}, reverse=True)
        self.category = ct_graph_in.category
        self.total_residuals = {}
        # TODO check if I truly need all three of msg new, msglog and msg. chech if i need both fullresidual and fullresidual new,

        for node in ct_graph_in.nodes(data=True):
            if node[1]["category"] == "factor":
                self.initial_beliefs[node[0]] = node[1]["InitialBelief"].Factor.array
            elif (
                    node[1]["category"] == "peptide"
                    or node[1]["category"] == ct_graph_in.category
            ):
                self.initial_beliefs[node[0]] = np.array(
                    [node[1]["InitialBelief_0"], node[1]["InitialBelief_1"]]
                )
            else:
                self.initial_beliefs[node[0]] = np.ones(
                    4
                )  # this entry will never be used as the convolution trees do not hold beliefs

        self.current_beliefs = self.initial_beliefs.copy()

        for node1, node2, data in ct_graph_in.edges(data=True):
            start_name, end_name = node1, node2

            if "MessageLength" in data:
                self.msg[(start_name, end_name)] = np.ones(data["MessageLength"])
            else:
                self.msg[(start_name, end_name)] = np.array([0.5, 0.5])

            self.msg[(end_name, start_name)] = self.msg[(start_name, end_name)]

            if "MessageLength" in data:
                self.msg_new[(start_name, end_name)] = np.ones(data["MessageLength"])
            else:
                self.msg_new[(start_name, end_name)] = np.array([0, 0])

            self.msg_new[(end_name, start_name)] = self.msg_new[(start_name, end_name)]

        self.msg_log = self.msg_new.copy()

    # variables (peptides,proteins,taxa)
    def get_incoming_message_variable(self, node, node_in):
        return self.msg[node, node_in]

    def compute_out_message_variable(self, node_out, node_in):
        incoming_messages = []
        node_belief = self.current_beliefs[node_out]
        for node_out_neighbor in self.graph.neighbors(node_out):
            if node_out_neighbor != node_in:
                incoming_messages.append(
                    self.get_incoming_message_variable(node_out_neighbor, node_out)
                )

        if not incoming_messages:
            # TODO Tanja: what exactly do you want to compare here?
            return node_belief if any(node_belief == self.initial_beliefs[node_out]) else self.msg[node_out, node_in]

        # need for logs to prevent underflow in very large multiplications
        incoming_messages = np.asarray(np.log(incoming_messages)).reshape(
            len(incoming_messages), 2
        )

        out_message_log = log_normalize(
            np.asarray(
                [
                    np.sum([np.log(node_belief[0]), np.sum(incoming_messages[:, 0])]),
                    np.sum([np.log(node_belief[1]), np.sum(incoming_messages[:, 1])]),
                ]
            )
        )

        if not np.all(out_message_log):
            out_message_log[out_message_log == 0] = 1e-30

        return out_message_log

    # factors (Conditional probability tables), handles different dimension of output/input variables
    def get_incoming_message_factor(self, node, node_in):
        return self.msg[node, node_in]

    def compute_out_message_factor(self, node_out, node_in):
        incoming_messages = []
        node_belief = self.current_beliefs[node_out]

        for node_out_neighbour in self.graph.neighbors(node_out):
            if node_out_neighbour != node_in:
                if [self.get_incoming_message_factor(node_out_neighbour, node_out)]:
                    # only the messages that have changed get multiplied into the current belief again
                    incoming_messages.append(
                        self.get_incoming_message_factor(node_out_neighbour, node_out)
                    )

        if self.graph.nodes[node_in]["category"] == "Convolution Tree":
            # handles empty & messages with only one value
            incoming_messages.append([1.0, 1.0])
            incoming_messages = np.asarray(incoming_messages).reshape(
                len(incoming_messages), 2
            )  # np.asarray(np.log(IncomingMessages)).reshape(len(IncomingMessages),2)
            out_messages = normalize(
                np.multiply(
                    node_belief,
                    [np.prod(incoming_messages[:, 0]), np.prod(incoming_messages[:, 1])],
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
                out_messages_log = log_normalize(np.add(log_belief, incoming_messages_log))
                if not np.all(out_messages_log):
                    out_messages_log[out_messages_log == 0] = 1e-30

                return [np.sum(out_messages_log[0, :]), np.sum(out_messages_log[1, :])]
            else:
                incoming_messages.append([1.0, 1.0])
                incoming_messages = np.asarray(incoming_messages).reshape(
                    len(incoming_messages), 2
                )
                out_messages = normalize(
                    np.multiply(
                        node_belief,
                        [
                            np.prod(incoming_messages[:, 0]),
                            np.prod(incoming_messages[:, 1]),
                        ],
                    )
                )
                return [np.sum(out_messages[0, :]), np.sum(out_messages[1, :])]

                # CTree, computes all out messages in one go

    def compute_out_messages_ct_tree(self, node):
        prot_prob_list = []
        old_prot_prob_list = []
        shared_likelihoods = np.ones(self.graph.nodes[node]["NumberOfParents"] + 1)
        peptides = []
        prot_list = []

        for node_in in self.graph.neighbors(node):
            if "CPD" not in str(node_in):
                if isinstance(self.msg[node_in, node], list):
                    prot_prob_list.append([self.msg[node_in, node]])
                else:
                    prot_prob_list.append([self.msg[node_in, node].tolist()])

                if isinstance(self.msg_log[node_in, node], list):
                    old_prot_prob_list.append([self.msg_log[node_in, node]])
                else:
                    old_prot_prob_list.append([self.msg_log[node_in, node].tolist()])

                prot_list.append(node_in)
            else:
                peptides.append(node_in)
                try:
                    shared_likelihoods = np.multiply(
                        avoid_underflow(shared_likelihoods), self.msg[node_in, node]
                    )
                except:
                    print(shared_likelihoods, self.msg[node_in, node])
                old_shared_likelihoods = np.multiply(
                    avoid_underflow(shared_likelihoods), self.msg_log[node_in, node]
                )

        # TODO Tanja: what are you trying to compare here? Do you want to check if all values in these lists
        # are different?
        if all(old_shared_likelihoods != shared_likelihoods) and any(
                [
                    (prot_prob_list[i][0]) != (old_prot_prob_list[i][0])
                    for i in range(len(prot_prob_list))
                ]
        ):
            # only update when the shared likelihoods or at least on of the protein messages has changed
            ct = ConvolutionTree(shared_likelihoods, prot_prob_list)

            for protein in range(len(prot_list)):
                self.msg_new[node, prot_list[protein]] = ct.MessageToVariable(protein)
                if not np.all(self.msg_new[node, prot_list[protein]]):
                    self.msg_new[node, prot_list[protein]][
                        self.msg_new[node, prot_list[protein]] == 0
                        ] = 1e-30

            for pep in peptides:
                self.msg_new[node, pep] = ct.MessageToSharedLikelihood()
                if not np.all(self.msg_new[node, pep]):
                    self.msg_new[node, pep][self.msg_new[node, pep] < 1e-30] = 1e-30

        else:
            for protein in range(len(prot_list)):
                self.msg_new[node, prot_list[protein]] = self.msg[node, prot_list[protein]]

            for pep in peptides:
                self.msg_new[node, pep] = self.msg[node, pep]

    # keeps track of which CTs have been update already in the current computeUpdate() loop
    def CTupdatecheck(self, CT):
        if CT in self.ListOfCTs:
            return False
        else:
            self.ListOfCTs.add(CT)
            return True

    # computes the residual between message new/message for a give edge(nodein/nodeOUT)
    def ComputeResidual(self, NodeIN, NodeOUT):
        Msg1 = self.msg_new[NodeIN, NodeOUT]
        Msg2 = self.msg[NodeIN, NodeOUT]
        print(Msg1, Msg2)
        if len(self.msg_new[NodeIN, NodeOUT]) != len(self.msg[NodeIN, NodeOUT]):
            Msg2 = [1] * len(self.msg_new[NodeIN, NodeOUT])
        return np.sum(np.abs(np.subtract(Msg1, Msg2)))

    def ComputeInfinityNormResidual(self, StartName, EndName):
        Msg1 = self.msg[StartName, EndName]
        Msg2 = self.msg_log[StartName, EndName]
        if len(self.msg_log[StartName, EndName]) != len(self.msg[StartName, EndName]):
            Msg2 = [1] * len(self.msg[StartName, EndName])
        pos = 0
        for i in Msg1:
            if i < 1e-30:
                Msg1[pos] = 1e-30
            pos += 1

        return np.max(np.abs(np.log(np.divide(Msg1, Msg2))))

        # approximate residual with zero look-ahead

    def ComputeZeroLookAheadResidual(self, StartName, EndName):
        NodeINneighbors = [nodes for nodes in self.graph.neighbors(StartName)]
        NodeINneighbors.remove(EndName)
        ApproximateResidual = sum(
            [self.full_residual[(neighbors, StartName)] for neighbors in NodeINneighbors]
        )
        return ApproximateResidual

    def ComputeTotalResiduals(self, StartName, EndName, CurrentResidual):
        for startNeighbors in self.graph.neighbors(StartName):
            if startNeighbors != EndName:
                self.total_residuals[
                    ((startNeighbors, StartName), (StartName, EndName))
                ] = 0

        for EndNeighbors in self.graph.neighbors(EndName):
            if EndNeighbors != StartName:
                self.total_residuals[((StartName, EndName), (EndName, EndNeighbors))] = (
                        self.total_residuals[((StartName, EndName), (EndName, EndNeighbors))]
                        + CurrentResidual
                )

    def compute_priority(self, start_name, end_name):
        self.priorities[(start_name, end_name)] = 0

        for end_neighbor in self.graph.neighbors(end_name):
            if end_neighbor != start_name:
                self.priorities[end_name, end_neighbor] = np.sum(
                    [
                        self.total_residuals[(sum_run, end_name), (end_name, end_neighbor)]
                        for sum_run in self.graph.neighbors(end_name)
                        if sum_run != end_neighbor
                    ]
                )

    # computes new message for a given edge (startname,endname) in the direction startname->endname
    def SingleEdgeDirectionUpdate(self, StartName, EndName):
        if (
                self.graph.nodes[StartName]["category"] == self.category
                or self.graph.nodes[StartName]["category"] == "peptide"
        ):
            self.msg_new[StartName, EndName] = self.compute_out_message_variable(
                StartName, EndName
            )

        if self.graph.nodes[StartName]["category"] == "Convolution Tree":
            CTCheck = self.CTupdatecheck(StartName)

            if CTCheck:
                self.compute_out_messages_ct_tree(StartName)

        if self.graph.nodes[StartName]["category"] == "factor":
            self.msg_new[StartName, EndName] = self.compute_out_message_factor(
                StartName, EndName
            )

    # compute updated messages for all edges
    def ComputeUpdate(self, localloops=False):
        self.ListOfCTs = set()  # keeps track of which CT has already been active

        if not isinstance(localloops, bool):
            raise TypeError("localloops needs to be boolean")

        if localloops and self.max_val:
            for EndName in self.graph.neighbors(self.max_val[1]):
                StartName = self.max_val[1]
                self.SingleEdgeDirectionUpdate(StartName, EndName)

        else:
            for edge in self.graph.edges():
                # update all edges
                StartName, EndName = edge[0], edge[1]
                self.SingleEdgeDirectionUpdate(StartName, EndName)

                StartName, EndName = edge[1], edge[0]
                self.SingleEdgeDirectionUpdate(StartName, EndName)

    def updateResidualMessage(self, Residual):
        """
        check which message Residual has the largest residual and updates that message in self.Msg
        :param Residual: dict, residuals of the last belief propagation iteration
        """

        self.max_val = max(Residual, key=Residual.get)
        self.msg[self.max_val] = self.msg_new[self.max_val]
        return Residual[self.max_val]

    def getPriorityMessage(self, PriorityVector):
        self.max_val = PriorityVector.top()
        return self.max_val

    def ZeroLookAheadLoopyLoop(self, max_loops, tolerance, local=False):
        """
        Run the zero-look-ahead belief propagation algorithm.
        :param max_loops: int, maximum number of iterations in case of non-convergence
        :param tolerance: float, tolerance for convergence check
        :param local: Bool, parameter passed to Computed Update function
        """

        max_residual = 100

        # first, do 5 loops where I update all messages
        print(f"Time spent in loop 0/{max_loops}: 0s", end="")
        for k in range(0, 5):
            start_t = time.time()
            self.ComputeUpdate()
            self.msg_log.update(self.msg)
            self.msg.update(self.msg_new)
            k += 1
            end_t = time.time()
            print(f"\rTime spent in loop {k}/{max_loops}: {(end_t - start_t):.3f}s", end="")

        # compute all residuals after 5 runs once (= initialize the residual/priorities vectors)
        for edge in self.graph.edges():
            # compute all residuals of the messages in this loop
            StartName, EndName = edge[1], edge[0]
            self.full_residual[
                (StartName, EndName)
            ] = self.ComputeInfinityNormResidual(StartName, EndName)

            # initialize the total residual to 0
            for End2 in self.graph.neighbors(EndName):
                self.total_residuals[((StartName, EndName), (EndName, End2))] = 0
                self.total_residuals[((End2, EndName), (EndName, StartName))] = 0

            StartName, EndName = edge[0], edge[1]
            self.full_residual[
                (StartName, EndName)
            ] = self.ComputeInfinityNormResidual(StartName, EndName)

            # initialize the total residual to 0
            for End2 in self.graph.neighbors(EndName):
                self.total_residuals[((StartName, EndName), (EndName, End2))] = 0
                self.total_residuals[((End2, EndName), (EndName, StartName))] = 0

        # set the priority vector once with copy of the previously calculated residuals
        self.priorities = pqdict(self.full_residual.copy(), reverse=True)

        k = 5

        print(f"\rTime spent in loop 0/{max_loops}: 0s -> residual max 0", end="")
        while k < max_loops and max_residual > tolerance:
            # actual zero-look-ahead-BP part
            start_t = time.time()
            priority_message = self.getPriorityMessage(self.priorities)
            max_residual = self.priorities[priority_message]
            self.SingleEdgeDirectionUpdate(priority_message[0], priority_message[1])
            priority_residual = self.ComputeInfinityNormResidual(
                priority_message[0], priority_message[1]
            )
            self.msg_log.update(self.msg)
            self.msg.update(self.msg_new)
            self.ComputeTotalResiduals(
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
        for variable in self.graph.nodes():
            if (
                    self.graph.nodes[variable]["category"] == self.category
                    or self.graph.nodes[variable]["category"] == "peptide"
            ):
                IncomingMessages = []

                for VariableNeighbors in self.graph.neighbors(variable):
                    IncomingMessages.append(
                        self.get_incoming_message_variable(VariableNeighbors, variable)
                    )

                # log to avoid overflow
                IncomingMessages = np.asarray(np.log(IncomingMessages)).reshape(
                    len(IncomingMessages), 2
                )
                LoggedVariableMarginal = log_normalize(
                    np.asarray(
                        [
                            np.sum(
                                [
                                    np.log(self.initial_beliefs[variable][0]),
                                    np.sum(IncomingMessages[:, 0]),
                                ]
                            ),
                            np.sum(
                                [
                                    np.log(self.initial_beliefs[variable][1]),
                                    np.sum(IncomingMessages[:, 1]),
                                ]
                            ),
                        ]
                    )
                )

                self.current_beliefs[variable] = LoggedVariableMarginal

    def DetectOscillations(self):
        pass


# calibration through message passing of all subgraphs in the List of factor graphs
def CalibrateAllSubgraphs(ListOfCTFactorGraphs, MaxIterations, Tolerance, local=False):
    """
    Performs bayesian inference through loopy belief propagation, returns dictionary {variable:posterior_probability}
    :param ListOfCTFactorGraphs: list, contains FactorGraph objects on which inference can be performed
    :param MaxIterations: int, max number of iterations in case of non-convergence
    :param Tolerance: float, error tolerance between messages for convergence criterion
    :param local: Bool, whether loops are calculated locally
    """

    if not isinstance(ListOfCTFactorGraphs, list):
        raise TypeError("ListOfFactorGraphs needs to be a list of graphs")
    if not isinstance(local, bool):
        raise TypeError("localloops needs to be boolean")

    ResultsList = []
    ResultsDict = {}
    NodeDict = {}

    for (idx, graph) in enumerate(ListOfCTFactorGraphs):
        if graph.number_of_nodes() > 2:
            NodeDict.update(dict(graph.nodes(data="category")))
            InitializedMessageObject = Messages(graph)
            InitializedMessageObject.ZeroLookAheadLoopyLoop(
                MaxIterations, Tolerance, local
            )
            ResultsList.append(InitializedMessageObject.current_beliefs)
            ResultsDict.update(InitializedMessageObject.current_beliefs)

    return ResultsList, ResultsDict, NodeDict


def SaveResultsToCsv(ResultsDict, NodeDict, NameString):
    """
    Save Loopy Belief Propagation results to .csv file
    :param ResultsDict: dict, {variable:posterior_probability}
    :param NodeDict: dict, dictionary of nodes that were in the factor graph and their attributes, to include the node category in the results
    :param NameString: str, csv output path
    """

    if not isinstance(NameString, str):
        raise TypeError("AddNameString needs to a string with Info on your run")
    if not isinstance(ResultsDict, dict):
        raise TypeError("Resultsdict must be dictionary")
    if not isinstance(NodeDict, dict):
        raise TypeError("Resultsdict must be dictionary")

    FullResultsDict = {
        key: [ResultsDict[key][1], NodeDict[key]] for key in ResultsDict.keys()
    }
    pd.DataFrame.from_dict(data=FullResultsDict, orient="index").to_csv(
        NameString, header=False
    )  # (datetime.now().strftime("%Y-%-m-%d-%H-%M-%S")+'-' +NameString + '-Results.csv', header = False)

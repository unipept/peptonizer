# implementation of belief propagation on a peptide-protein graph
# __________________________________________________________________________________________
import numpy as np
import networkx as nx
import math
import pandas as pd

from scipy.signal import fftconvolve
from factor_graph_generation import *
from scipy.special import logsumexp

import time


# set numpy to raise errors
# np.seterr(all='raise')


def normalize(array):
    return array / np.sum(array)


# normalization of log propabilities
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
        resultcheck = fftconvolve(otherJointVector[::-1], self.likelihoodBelow)
        result = fftconvolve(otherJointVector[::-1], self.likelihoodBelow)[
            startingPoint : startingPoint + answerSize
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
        for i in range(0, 2**self.logLength - len(proteins)):
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
        return self.LastNode.jointAbove[0][0 : (self.nProteins + 1)]


class Messages:
    """
    Class holding all messages and beliefs.
    Functions execute loopy residual belief propagation with zero look ahead
    """

    # class that holds the messages of itereation t and iteration t+1 as dictionaries

    def __init__(self, CTGraphIn):
        self.Msg = {}
        self.MsgNew = {}
        self.MsgLog = {}
        self.Graph = CTGraphIn
        self.MaxVal = None
        self.FullResidual = {}
        self.FullResidualNew = {}
        self.InitialBeliefs = {}
        self.CurrentBeliefs = {}
        self.CurrentBeliefsNew = {}
        self.queue = {}
        self.priorities = {}
        self.category = CTGraphIn.category
        self.TotalResiduals = {}
        # TODO check if I truly need all three of msg new, msglog and msg. chech if i need both fullresidual and fullresidual new,

        for node in CTGraphIn.nodes(data=True):
            if node[1]["category"] == "factor":
                self.InitialBeliefs[node[0]] = node[1]["InitialBelief"].Factor.array
            elif (
                node[1]["category"] == "peptide"
                or node[1]["category"] == CTGraphIn.category
            ):
                self.InitialBeliefs[node[0]] = np.array(
                    [node[1]["InitialBelief_0"], node[1]["InitialBelief_1"]]
                )
            else:
                self.InitialBeliefs[node[0]] = np.ones(
                    4
                )  # this entry will never be used as the convolution trees do not hold beliefs

        self.CurrentBeliefs = self.InitialBeliefs.copy()
        self.CurrentBeliefsNew = self.InitialBeliefs.copy()

        for node1, node2, data in CTGraphIn.edges(data=True):
            StartName, EndName = node1, node2

            if "MessageLength" in data:
                self.Msg[(StartName, EndName)] = np.ones(data["MessageLength"])
            else:
                self.Msg[(StartName, EndName)] = np.array([0.5, 0.5])

            self.Msg[(EndName, StartName)] = self.Msg[(StartName, EndName)]

            if "MessageLength" in data:
                self.MsgNew[(StartName, EndName)] = np.ones(data["MessageLength"])
            else:
                self.MsgNew[(StartName, EndName)] = np.array([0, 0])

            self.MsgNew[(EndName, StartName)] = self.MsgNew[(StartName, EndName)]

        self.MsgLog = self.MsgNew.copy()

    # variables (peptides,proteins,taxa)
    def GetIncomingMessageVariable(self, Node, NodeIN):
        returnedMessage = self.Msg[Node, NodeIN]
        return returnedMessage

    def ComputeOutMessageVariable(self, NodeOUT, NodeIN):
        IncomingMessages = []
        NodeBelief = self.CurrentBeliefs[NodeOUT]
        for NodeOUTneighbors in self.Graph.neighbors(NodeOUT):
            if NodeOUTneighbors != NodeIN:
                IncomingMessages.append(
                    self.GetIncomingMessageVariable(NodeOUTneighbors, NodeOUT)
                )

        if not IncomingMessages:
            check = any(NodeBelief == self.InitialBeliefs[NodeOUT])
            if any(NodeBelief == self.InitialBeliefs[NodeOUT]):
                return NodeBelief
            else:
                return self.Msg[NodeOUT, NodeIN]

        else:
            # need for logs to prevent underflow in very large multiplications
            IncomingMessages = np.asarray(np.log(IncomingMessages)).reshape(
                len(IncomingMessages), 2
            )
            OutMessageLog = log_normalize(
                np.asarray(
                    [
                        np.sum([np.log(NodeBelief[0]), np.sum(IncomingMessages[:, 0])]),
                        np.sum([np.log(NodeBelief[1]), np.sum(IncomingMessages[:, 1])]),
                    ]
                )
            )
            if np.isnan(np.sum(OutMessageLog)) == True:
                stoppoint = 3
            if not np.all(OutMessageLog):
                OutMessageLog[OutMessageLog == 0] = 1e-30
                stoppoint = 5
            return OutMessageLog

    # factors (Conditional probability tables), handles different dimension of output/input variables
    def GetIncomingMessageFactor(self, Node, NodeIN):
        check1 = self.Msg[Node, NodeIN]
        check2 = self.MsgLog[Node, NodeIN]
        returnedMessage = self.Msg[Node, NodeIN]
        return returnedMessage

    def ComputeOutMessageFactor(self, NodeOUT, NodeIN):
        IncomingMessages = []
        NodeBelief = self.CurrentBeliefs[NodeOUT]

        for NodeOUTneighbors in self.Graph.neighbors(NodeOUT):
            if NodeOUTneighbors != NodeIN:
                if [self.GetIncomingMessageFactor(NodeOUTneighbors, NodeOUT)]:
                    # only the messages that have changed get multiplied into the current belief again
                    IncomingMessages.append(
                        self.GetIncomingMessageFactor(NodeOUTneighbors, NodeOUT)
                    )

        if self.Graph.nodes[NodeIN]["category"] == "Convolution Tree":
            # handles empty & messages with only one value
            IncomingMessages.append([1.0, 1.0])
            IncomingMessages = np.asarray(IncomingMessages).reshape(
                len(IncomingMessages), 2
            )  # np.asarray(np.log(IncomingMessages)).reshape(len(IncomingMessages),2)
            OutMessages = normalize(
                np.multiply(
                    NodeBelief,
                    [np.prod(IncomingMessages[:, 0]), np.prod(IncomingMessages[:, 1])],
                )
            )  # lognormalize(np.add(np.log(NodeBelief),[np.sum(IncomingMessages[:,0]),np.sum(IncomingMessages[:,1])]))#
            # self.CurrentBeliefsNew[NodeOUT] = OutMessages

            return np.add(OutMessages[:, 0], OutMessages[:, 1])
        else:
            if np.asarray(IncomingMessages[0]).shape[0] > 2:
                IncomingMessagesShaped = np.asarray(IncomingMessages).reshape(
                    IncomingMessages[0].shape[0], 1
                )
                # OutMessages = normalize(NodeBelief*IncomingMessagesShaped)

                IncomingMessagesLog = np.log(
                    np.asarray(IncomingMessages).reshape(
                        IncomingMessages[0].shape[0], 1
                    )
                )
                # catch warning for log(0)

                LogBelief = np.log(NodeBelief)
                OutMessagesLog = log_normalize(np.add(LogBelief, IncomingMessagesLog))
                # if not np.all([np.sum(OutMessages[0,:]),np.sum(OutMessages[1,:])]):
                #    stoopoint = 4
                if not np.all(OutMessagesLog):
                    OutMessagesLog[OutMessagesLog == 0] = 1e-30
                    stooptpoint = 5

                return [np.sum(OutMessagesLog[0, :]), np.sum(OutMessagesLog[1, :])]
            else:
                IncomingMessages.append([1.0, 1.0])
                IncomingMessages = np.asarray(IncomingMessages).reshape(
                    len(IncomingMessages), 2
                )
                OutMessages = normalize(
                    np.multiply(
                        NodeBelief,
                        [
                            np.prod(IncomingMessages[:, 0]),
                            np.prod(IncomingMessages[:, 1]),
                        ],
                    )
                )
                return [np.sum(OutMessages[0, :]), np.sum(OutMessages[1, :])]

                # CTree, computes all out messages in one go

    def ComputeOutMessagesCTtree(self, Node):
        ProtProbList = []
        OldProtProbList = []
        sharedLikelihoods = np.ones(self.Graph.nodes[Node]["NumberOfParents"] + 1)
        peptides = []
        ProtList = []
        self.CurrentBeliefsNew[Node] = np.ones(4)

        for NodesIN in self.Graph.neighbors(Node):
            if "CPD" not in str(NodesIN):
                if type(self.Msg[NodesIN, Node]) == list:
                    ProtProbList.append([self.Msg[NodesIN, Node]])
                else:
                    ProtProbList.append([self.Msg[NodesIN, Node].tolist()])

                if type(self.MsgLog[NodesIN, Node]) == list:
                    OldProtProbList.append([self.MsgLog[NodesIN, Node]])
                else:
                    OldProtProbList.append([self.MsgLog[NodesIN, Node].tolist()])

                ProtList.append(NodesIN)
            else:
                peptides.append(NodesIN)
                try:
                    sharedLikelihoods = np.multiply(
                        avoid_underflow(sharedLikelihoods), self.Msg[NodesIN, Node]
                    )
                except:
                    print(sharedLikelihoods, self.Msg[NodesIN, Node])
                OldSharedLikelihoods = np.multiply(
                    avoid_underflow(sharedLikelihoods), self.MsgLog[NodesIN, Node]
                )

        if all(OldSharedLikelihoods != sharedLikelihoods) and any(
            [
                (ProtProbList[i][0]) != (OldProtProbList[i][0])
                for i in range(len(ProtProbList))
            ]
        ):
            # only update when the shared likelihhods or at least on of the protein messages has changed
            CT = ConvolutionTree(sharedLikelihoods, ProtProbList)

            for protein in range(len(ProtList)):
                self.MsgNew[Node, ProtList[protein]] = CT.MessageToVariable(protein)
                if not np.all(self.MsgNew[Node, ProtList[protein]]):
                    self.MsgNew[Node, ProtList[protein]][
                        self.MsgNew[Node, ProtList[protein]] == 0
                    ] = 1e-30
                    stop = 3

            for pep in peptides:
                self.MsgNew[Node, pep] = CT.MessageToSharedLikelihood()
                if not np.all(self.MsgNew[Node, pep]):
                    self.MsgNew[Node, pep][self.MsgNew[Node, pep] < 1e-30] = 1e-30
                    stopp = 5

        else:
            for protein in range(len(ProtList)):
                self.MsgNew[Node, ProtList[protein]] = self.Msg[Node, ProtList[protein]]

            for pep in peptides:
                self.MsgNew[Node, pep] = self.Msg[Node, pep]

    # keeps track of which CTs have been update already in the current computeUpdate() loop
    def CTupdatecheck(self, CT):
        if CT in self.ListOfCTs:
            return False
        else:
            self.ListOfCTs.append(CT)
            return True

    # computes the residual between message new/message for a give edge(nodein/nodeOUT)
    def ComputeResidual(self, NodeIN, NodeOUT):
        Msg1 = self.MsgNew[NodeIN, NodeOUT]
        Msg2 = self.Msg[NodeIN, NodeOUT]
        print(Msg1, Msg2)
        if len(self.MsgNew[NodeIN, NodeOUT]) != len(self.Msg[NodeIN, NodeOUT]):
            Msg2 = [1] * len(self.MsgNew[NodeIN, NodeOUT])
        return np.sum(np.abs(np.subtract(Msg1, Msg2)))

    def ComputeInfinityNormResidual(self, StartName, EndName):
        Msg1 = self.Msg[StartName, EndName]
        Msg2 = self.MsgLog[StartName, EndName]
        if len(self.MsgLog[StartName, EndName]) != len(self.Msg[StartName, EndName]):
            Msg2 = [1] * len(self.Msg[StartName, EndName])
        check = np.max(np.abs(np.log(np.divide(Msg1, Msg2))))
        pos = 0
        for i in Msg1:
            if i < 1e-30:
                Msg1[pos] == 1e-30
            pos += 1

        return np.max(np.abs(np.log(np.divide(Msg1, Msg2))))

        # approximate residual with zero look-ahead

    def ComputeZeroLookAheadResidual(self, StartName, EndName):
        NodeINneighbors = [nodes for nodes in self.Graph.neighbors(StartName)]
        NodeINneighbors.remove(EndName)
        ApproximateResidual = sum(
            [self.FullResidual[(neighbors, StartName)] for neighbors in NodeINneighbors]
        )
        return ApproximateResidual

    def ComputeTotalResiduals(self, StartName, EndName, CurrentResidual):
        for startNeighbors in self.Graph.neighbors(StartName):
            if startNeighbors != EndName:
                self.TotalResiduals[
                    ((startNeighbors, StartName), (StartName, EndName))
                ] = 0

        for EndNeighbors in self.Graph.neighbors(EndName):
            if EndNeighbors != StartName:
                check = self.TotalResiduals[
                    ((StartName, EndName), (EndName, EndNeighbors))
                ]
                self.TotalResiduals[((StartName, EndName), (EndName, EndNeighbors))] = (
                    self.TotalResiduals[((StartName, EndName), (EndName, EndNeighbors))]
                    + CurrentResidual
                )

    def ComputePriority(self, StartName, EndName):
        self.priorities[(StartName, EndName)] = 0

        for EndNeighbors in self.Graph.neighbors(EndName):
            if EndNeighbors != StartName:
                check = np.sum(
                    [
                        self.TotalResiduals[(SumRun, EndName), (EndName, EndNeighbors)]
                        for SumRun in self.Graph.neighbors(EndName)
                        if SumRun != EndNeighbors
                    ]
                )
                self.priorities[EndName, EndNeighbors] = np.sum(
                    [
                        self.TotalResiduals[(SumRun, EndName), (EndName, EndNeighbors)]
                        for SumRun in self.Graph.neighbors(EndName)
                        if SumRun != EndNeighbors
                    ]
                )

        # for startNeighbors in self.Graph.neighbors(StartName):
        #    if startNeighbors != EndName:
        #        self.priorities[startNeighbors,EndName] = np.sum([self.TotalResiduals[(SumRun,EndName),(EndName,EndNeighbors)] for SumRun in self.Graph.neighbors(EndName) if SumRun != EndNeighbors])

    # computes new message for a given edge (startname,endname) in the direction startname->endname
    def SingleEdgeDirectionUpdate(self, StartName, EndName):
        if (
            self.Graph.nodes[StartName]["category"] == self.category
            or self.Graph.nodes[StartName]["category"] == "peptide"
        ):
            self.MsgNew[StartName, EndName] = self.ComputeOutMessageVariable(
                StartName, EndName
            )

        if self.Graph.nodes[StartName]["category"] == "Convolution Tree":
            CTCheck = self.CTupdatecheck(StartName)

            if CTCheck:
                self.ComputeOutMessagesCTtree(StartName)

        if self.Graph.nodes[StartName]["category"] == "factor":
            self.MsgNew[StartName, EndName] = self.ComputeOutMessageFactor(
                StartName, EndName
            )

    # compute updated messages for all edges
    def ComputeUpdate(self, localloops=False):
        self.ListOfCTs = []  # keeps track of which CT has already been active

        if not isinstance(localloops, bool):
            raise TypeError("localloops needs to be boolean")

        if localloops and self.MaxVal:
            for EndName in self.Graph.neighbors(self.MaxVal[1]):
                StartName = self.MaxVal[1]
                self.SingleEdgeDirectionUpdate(StartName, EndName)

        else:
            for edge in self.Graph.edges():
                # update all edges
                StartName, EndName = edge[0], edge[1]
                self.SingleEdgeDirectionUpdate(StartName, EndName)

                StartName, EndName = edge[1], edge[0]
                self.SingleEdgeDirectionUpdate(StartName, EndName)

        # for edge in self.Graph.edges():
        #    #compute all residuals of the messages in this loop
        #    StartName, EndName = edge[1], edge[0]
        #    self.FullResidual[(StartName, EndName)] = self.ComputeResidual(StartName, EndName)

        #    StartName, EndName = edge[0], edge[1]
        #    self.FullResidual[(StartName, EndName)] = self.ComputeResidual(StartName, EndName)

    def updateResidualMessage(self, Residual):
        """
        check which message Residual has the largest residual and updates that message in self.Msg
        :param Residual: dict, residuals of the last belief propagation iteration
        """

        self.MaxVal = max(Residual, key=Residual.get)
        self.Msg[self.MaxVal] = self.MsgNew[self.MaxVal]
        return Residual[self.MaxVal]

    def getPriorityMessage(self, PriorityVector):
        self.Maxval = max(PriorityVector, key=PriorityVector.get)
        return max(PriorityVector, key=PriorityVector.get)

    def ZeroLookAheadLoopyLoop(self, max_loops, tolerance, local=False):
        """
        Run the zero-look-ahead belief propagation algorithm.
        :param max_loops: int, maximum number of iterations in case of non-convergence
        :param tolerance: float, tolerance for convergence check
        :param local: Bool, parameter passed to Computed Update function
        """

        max_residual = 100

        # first, do 5 loops where I update all messages
        print("time per loop 0 0", end="")
        for k in range(0, 5):
            start_t = time.time()
            self.ComputeUpdate()
            self.MsgLog.update(self.Msg)
            self.Msg.update(self.MsgNew)
            k += 1
            end_t = time.time()
            print(f"\rtime per loop {k:.3f}s {end_t - start_t}", end="")

        # compute all residuals after 5 runs once (= initialize the residual/priorities vectors)
        for edge in self.Graph.edges():
            # compute all residuals of the messages in this loop
            StartName, EndName = edge[1], edge[0]
            self.FullResidual[
                (StartName, EndName)
            ] = self.ComputeInfinityNormResidual(StartName, EndName)

            # initialize the total residual to 0
            for End2 in self.Graph.neighbors(EndName):
                self.TotalResiduals[((StartName, EndName), (EndName, End2))] = 0
                self.TotalResiduals[((End2, EndName), (EndName, StartName))] = 0

            StartName, EndName = edge[0], edge[1]
            self.FullResidual[
                (StartName, EndName)
            ] = self.ComputeInfinityNormResidual(StartName, EndName)

            # set the priority vector once with copy of the previously calculated residuals
            self.priorities = self.FullResidual.copy()

            # initialize the total residual to 0
            for End2 in self.Graph.neighbors(EndName):
                self.TotalResiduals[((StartName, EndName), (EndName, End2))] = 0
                self.TotalResiduals[((End2, EndName), (EndName, StartName))] = 0

        k = 5

        print("time per loop 0 0 -> residual max 0")
        while k < max_loops and max_residual > tolerance:
            # actual zero-look-ahead-BP part
            start_t = time.time()
            # print(self.getPriorityMessage(self.priorities))
            priority_message = self.getPriorityMessage(self.priorities)
            max_residual = self.priorities[priority_message]
            self.SingleEdgeDirectionUpdate(priority_message[0], priority_message[1])
            priority_residual = self.ComputeInfinityNormResidual(
                priority_message[0], priority_message[1]
            )
            self.MsgLog.update(self.Msg)
            self.Msg.update(self.MsgNew)
            self.ComputeTotalResiduals(
                priority_message[0], priority_message[1], priority_residual
            )
            self.ComputePriority(priority_message[0], priority_message[1])

            end_t = time.time()

            # Only update the time per loop every 5 iterations
            if k % 5 == 0:
                print(f"\rtime per loop {k} {(end_t - start_t):.3f} -> residual max {max_residual:.3f}", end="")

            k += 1

        # marginalize once the model has converged
        for variable in self.Graph.nodes():
            if (
                self.Graph.nodes[variable]["category"] == self.category
                or self.Graph.nodes[variable]["category"] == "peptide"
            ):
                IncomingMessages = []

                for VariableNeighbors in self.Graph.neighbors(variable):
                    IncomingMessages.append(
                        self.GetIncomingMessageVariable(VariableNeighbors, variable)
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
                                    np.log(self.InitialBeliefs[variable][0]),
                                    np.sum(IncomingMessages[:, 0]),
                                ]
                            ),
                            np.sum(
                                [
                                    np.log(self.InitialBeliefs[variable][1]),
                                    np.sum(IncomingMessages[:, 1]),
                                ]
                            ),
                        ]
                    )
                )

                self.CurrentBeliefs[variable] = LoggedVariableMarginal

    def DetectOscillations(self):
        pass


# calibration through message passing of all subgraphs in the List of factor graphs
def CalibrateAllSubgraphs(ListOfCTFactorGraphs, MaxIterations, Tolerance, local=False):
    """
    Performs bayesian inferencethrough loopy belief propgation, returns dictionary {variable:posterior_probability}
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

    for Graph in ListOfCTFactorGraphs:
        if Graph.number_of_nodes() > 2:
            NodeDict.update(dict(Graph.nodes(data="category")))
            InitializedMessageObject = Messages(Graph)
            InitializedMessageObject.ZeroLookAheadLoopyLoop(
                MaxIterations, Tolerance, local
            )
            ResultsList.append(InitializedMessageObject.CurrentBeliefs)
            ResultsDict.update(InitializedMessageObject.CurrentBeliefs)

    return ResultsList, ResultsDict, NodeDict


# save the resulstsdictionary from CalibrateAllSubgraphs to a csv file
def SaveResultsToCsv(ResultsDict, NodeDict, NameString):
    """
    Save Loopy Belief Propagation results to .csv file
    :param ResultsDict: dict, {variable:posterior_probability}
    :param NodeDict: dict, dictionary of nodes that were in the factor graph and their attributes, to include the node category in the results
    :param NameString: str, csv ouptut path
    """

    if not isinstance(NameString, str):
        raise TypeError("AddNameString needs to a string with Info on your run")
    if not isinstance(ResultsDict, dict):
        raise TypeError("Resultsdict must be diciontary")
    if not isinstance(NodeDict, dict):
        raise TypeError("Resultsdict must be dictionary")

    FullResultsDict = {
        key: [ResultsDict[key][1], NodeDict[key]] for key in ResultsDict.keys()
    }
    pd.DataFrame.from_dict(data=FullResultsDict, orient="index").to_csv(
        NameString, header=False
    )  # (datetime.now().strftime("%Y-%-m-%d-%H-%M-%S")+'-' +NameString + '-Results.csv', header = False)
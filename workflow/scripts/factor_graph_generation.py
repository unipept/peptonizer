import numpy as np
import networkx as nx
from collections import namedtuple
import pandas as pd

def normalize(array):
    return array / np.sum(array)


def avoid_underflow(array):
    array[array < 1e-30] = 1e-30
    return array


class ProteinPeptideGraph(nx.Graph):
    # Class for storing the protein-peptide graph with scores and priors but no factors, e.g for visual representation.
    pass


class TaxonGraph(nx.Graph):
    """
    Class with functions to construct a peptide-taxon graph using Entrez/NCBI mapping.
    """

    def __init__(self):
        nx.Graph.__init__(self)
        self.TaxidList = []

    def CreateFromUnipeptResponseCSV(self, CSVpath):
        UnipeptResponse = pd.read_csv(CSVpath)
        # drop rows that have an entry in Highertaxa that appears only once
        counts = UnipeptResponse["HigherTaxa"].value_counts()
        UnipeptResponse = UnipeptResponse[
            UnipeptResponse["HigherTaxa"].isin(counts[counts > 1].index)
        ]
        newGraph = nx.from_pandas_edgelist(UnipeptResponse, "sequence", "HigherTaxa")
        PeptideAttributes = UnipeptResponse.apply(
            lambda row: (
                row["sequence"],
                {
                    "InitialBelief_0": row["score"],
                    "InitialBelief_1": 1 - row["score"],
                    "category": "peptide",
                },
            ),
            axis=1,
        )
        TaxaAttributes = UnipeptResponse.apply(
            lambda row: (row["HigherTaxa"], {"category": "taxon"}), axis=1
        )
        IntermediateGraph = nx.Graph()
        IntermediateGraph.add_edges_from(newGraph.edges)
        IntermediateGraph.add_nodes_from(PeptideAttributes)
        IntermediateGraph.add_nodes_from(TaxaAttributes)

        # cluster the resulting graph with the louvain algorithm
        communities = nx.community.louvain_communities(IntermediateGraph)
        # separate the graph into its communities and enter into same graph object
        for i, community in enumerate(communities):
            Subgraph = IntermediateGraph.subgraph(community)
            self.add_edges_from(Subgraph.edges)

        self.add_nodes_from(PeptideAttributes)
        self.add_nodes_from(TaxaAttributes)


class Factor:
    # represents noisy OR cpds, has dimension n(parensports)xn(peptide states(=2))
    def __init__(self, CPDarray, VariableArray):
        if isinstance(VariableArray, str):
            raise TypeError(
                "VariableArray: Expected type list or array like, got string"
            )

        Factor = namedtuple("Factor", ["array", "arrayLabels"])
        self.Factor = Factor(CPDarray, VariableArray)


class Variable:
    # has dimension of petide states ergo 2
    def __init__(self, ProbabilityArray, VariableArray):
        self.Variable = namedtuple("Variable", ["array", "arrayLabels"])
        Variable.array = ProbabilityArray
        Variable.arrayLabels = VariableArray


# the variable and factor types might be unecessary as i do not need a lot of flexibility in the input

# TODO implement checking& error raising if input aren't np.array/ list or array like


class FactorGraph(nx.Graph):
    def __init__(self):
        super().__init__()

    def ConstructFromProteinPeptideGraph(self, ProteinPeptideGraph):
        """'
        Takes a graph of proteins to peptides as input and adds the noisy-OR factors
        """
        nodelist = list(ProteinPeptideGraph.nodes(data=True))
        self.add_nodes_from(nodelist)
        for node in nodelist:
            # create noisy OR cpd per peptide
            if node[1]["category"] == "peptide":
                degree = ProteinPeptideGraph.degree(node[0])
                neighbors = list(ProteinPeptideGraph.neighbors(node[0]))
                self.add_node(node[0] + " CPD", category="factor", ParentNumber=degree)
                self.add_edges_from([(node[0] + " CPD", x) for x in neighbors])
                self.add_edge(node[0] + " CPD", node[0])

        return [self]

    def ConstructFromTaxonGraph(self, TaxonPeptideGraph):
        """'
        Takes a graph of Taxa to peptides in networkx form as input and adds the factor nodes
        """
        nodelist = list(TaxonPeptideGraph.nodes(data=True))
        self.add_nodes_from(nodelist)
        for node in nodelist:
            # create noisy OR cpd per peptide
            if node[1]["category"] == "peptide":
                degree = TaxonPeptideGraph.degree(node[0])
                neighbors = list(TaxonPeptideGraph.neighbors(node[0]))
                self.add_node(node[0] + " CPD", category="factor", ParentNumber=degree)
                self.add_edges_from([(node[0] + " CPD", x) for x in neighbors])
                self.add_edge(node[0] + " CPD", node[0])


# separate the connected components in the subgraph
def SeparateSubgraphs(graphIN, NodesToKeep):
    """
    separations of subgraphs (create news graphs for each subgraph)

    """
    newG = CTFactorGraph(nx.Graph())

    # add nodes to keep
    newG.add_nodes_from((n, graphIN.nodes[n]) for n in NodesToKeep)
    # add edges if they are present in original graph
    newG.add_edges_from(
        (n, nbr, d)
        for n, nbrs in graphIN.adj.items()
        if n in NodesToKeep
        for nbr, d in nbrs.items()
        if nbr in NodesToKeep
    )

    return newG  # ListOfFactorGraphs


class CTFactorGraph(FactorGraph):
    """'
    This class is a networkx graph representing the full graphical model with all variables, CTrees, and Noisy-OR factors
    """

    def __init__(self, GraphIn, GraphType="Taxons"):
        """
        takes either a graph or a path to a graphML file as input

        """
        super().__init__()

        GraphTypes = ["Proteins", "Taxons"]
        if GraphType not in GraphTypes:
            raise ValueError("Invalid Graphtype. Expected one of: %s" % GraphTypes)

        if GraphType == "Taxons":
            self.category = "taxon"
        elif GraphType == "Protein":
            self.category = "protein"

        if isinstance(GraphIn, str):
            GraphIn = nx.read_graphml(GraphIn)

        # need these to create a new instance of a CT fractorgraph and not overwrite the previous graph
        self.add_edges_from(GraphIn.edges(data=True))
        self.add_nodes_from(GraphIn.nodes(data=True))

    def AddCTNodes(self):
        """
        When creating the CTGraph and not just reading from a previously saved graph format, use this command to add the CT nodes
        """
        # create the convolution tree nodes and connect them in the graph
        ListOfEdgeAddList = []
        ListOfEdgeRemoveList = []
        ListOfProtLists = []
        ListOfCTs = []
        ListOfFactors = []
        for node in self.nodes(data=True):
            # go through all factors with degree>2 and get their protein lists, then generate their conv. trees
            if node[1]["category"] == "factor" and self.degree[node[0]] > 2:
                ProtList = []

                for neighbor in self.neighbors(node[0]):
                    neighbornode = self.nodes[neighbor]
                    if neighbornode["category"] == self.category:
                        ProtList.append(neighbor)

                ListOfCTs.append([1])
                ListOfFactors.append(node[0])
                ListOfProtLists.append(ProtList)
                ListOfEdgeAddList.append(
                    [("CTree " + " ".join(str(ProtList)), x) for x in ProtList]
                )
                ListOfEdgeRemoveList.append([(node[0], x) for x in ProtList])

        # Fill all info into graph structure, should probably do this inside the loop before, so that i can initialize the messages
        for i in range(len(ListOfCTs)):
            self.add_node(
                "CTree " + " ".join(str(ListOfProtLists[i])),
                ConvolutionTree=str(ListOfCTs[i][0]),
                category="Convolution Tree",
                NumberOfParents=len(ListOfProtLists[i]),
            )
            self.add_edge(
                "CTree " + " ".join(str(ListOfProtLists[i])),
                ListOfFactors[i],
                MessageLength=len(ListOfProtLists[i]) + 1,
            )
            self.add_edges_from(ListOfEdgeAddList[i])
            self.remove_edges_from(ListOfEdgeRemoveList[i])

    def SaveToGraphML(self, Filename):
        nx.write_graphml(self, Filename)

    def ComputeNetworkAttributes(self):
        """
        Computes nodes attributes using builtin networkx functions
        Returns degree centrality, closeness centrality, betweenness centrality and eigen centrality
        """
        DegreeCentrality = dict(
            sorted(nx.degree_centrality(self).items(), key=lambda item: item[1])
        )
        Closenesscentrality = dict(
            sorted(nx.closeness_centrality(self).items(), key=lambda item: item[1])
        )
        BetweennessCentrality = dict(
            nx.betweenness_centrality(self).items(), key=lambda item: item[1]
        )
        Eigencentrality = dict(
            nx.eigenvector_centrality(self).items(), key=lambda item: item[1]
        )

        return (
            DegreeCentrality,
            Closenesscentrality,
            BetweennessCentrality,
            Eigencentrality,
        )

    def FillInFactors(self, alpha, beta, regularized):
        """fills in the noisy or Factors according to detection and error probabilities given"""

        for node in self.nodes(data=True):
            # create noisy OR cpd per peptide
            if node[1]["category"] == "factor":
                # add noisyOR factors
                degree = node[1]["ParentNumber"]
                # pre-define the CPD array and fill it with the noisyOR values
                cpdArray = np.full([2, degree + 1], 1 - alpha)
                cpdArray_regularized = np.full([2, degree + 1], 1 - alpha)
                ExponentArray = np.arange(0, degree + 1)
                divideArray = np.concatenate(
                    (np.asarray([1]), np.arange(1, degree + 1))
                )
                # regularize cpd priors to penalize higher number of parents
                # log domain to avoid underflow
                cpdArray[0, :] = np.power(cpdArray[0, :], ExponentArray) * (1 - beta)
                cpdArray_regularized[0, :] = np.divide(
                    np.power(cpdArray[0, :], ExponentArray) * (1 - beta), divideArray
                )
                # check0 = cpdArray_regularized[0,:]
                # check1 = cpdArray[0,:]
                cpdArray[1, :] = np.add(-cpdArray[0, :], 1)
                cpdArray_regularized[1, :] = np.add(-cpdArray_regularized[0, :], 1)
                cpdArray = np.transpose(normalize(cpdArray))
                cpdArray_regularized = avoid_underflow(
                    np.transpose(normalize(cpdArray_regularized))
                )
                if regularized == True:
                    FactorToAdd = Factor(
                        cpdArray_regularized,
                        ["placeholder", [node[0] + "0", node[0] + "1"]],
                    )
                else:
                    FactorToAdd = Factor(
                        cpdArray, ["placeholder", [node[0] + "0", node[0] + "1"]]
                    )

                # add factor & its edges to network as an extra node
                nx.set_node_attributes(self, {node[0]: FactorToAdd}, "InitialBelief")

    def FillInPriors(self, prior):
        """fills in the taxon priors according to the given prior"""

        for node in self.nodes(data=True):
            # create noisy OR cpd per peptide
            if node[1]["category"] == "taxon":
                nx.set_node_attributes(
                    self,
                    {node[0]: {"InitialBelief_0": 1 - prior, "InitialBelief_1": prior}},
                )


def GenerateCTFactorGraphs(ListOfFactorGraphs, GraphType="Taxons"):
    if type(ListOfFactorGraphs) is not list:
        ListOfFactorGraphs = [ListOfFactorGraphs]
    for Graph in ListOfFactorGraphs:
        CTFactorgraph = CTFactorGraph(
            Graph, GraphType
        )  # ListOfCTFactorGraphs.append(CTFactorGraph(Graph,GraphType))
    return CTFactorgraph

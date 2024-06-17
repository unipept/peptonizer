import numpy as np
import math
import array_utils

from scipy.signal import fftconvolve


# implementation of the convolution tree according to serang
# (see: https://bitbucket.org/orserang/convolutiontree/src/master/SplicoformSolver.py)
# not written by me!!
class CTNode:
    def __init__(self, joint_above):
        # normalize for greater precision
        self.joint_above = array_utils.normalize(joint_above)

        self.left_parent = None
        self.right_parent = None

        self.likelihood_below = None

    # passing msges down: adding variables
    @classmethod
    def create_count_node(cls, lhs, rhs):
        # create a count node with joint prob for two parents above(vs the init if we have no parents)
        joint_above = fftconvolve(lhs.joint_above, rhs.joint_above)
        result = cls(joint_above)

        result.left_parent = lhs
        result.right_parent = rhs

        return result

    # passing messages up : subtracting variables
    def message_up(self, answer_size, other_joint_vector):
        starting_point = len(other_joint_vector) - 1
        result = fftconvolve(
            other_joint_vector[::-1],
            self.likelihood_below
        )[starting_point: starting_point + answer_size]
        return array_utils.normalize(result)

    def message_up_left(self):
        return self.message_up(
            len(self.left_parent.joint_above[0]), self.right_parent.joint_above[0]
        )

    def message_up_right(self):
        return self.message_up(
            len(self.right_parent.joint_above[0]), self.left_parent.joint_above[0]
        )

    # once all messages are received
    def posterior(self):
        return array_utils.normalize(self.joint_above * self.likelihood_below)

    def messages_up(self):
        return self.likelihood_below


class ConvolutionTree:
    def __init__(self, n_to_shared_likelihoods, proteins):
        self.n_to_shared_likelihoods = n_to_shared_likelihoods
        self.log_length = int(math.ceil(np.log2(float(len(proteins)))))  # length we need
        self.all_layers = []
        self.build_first_layer(proteins)
        self.build_remaining_layers()
        self.propagate_backward()
        self.n_proteins = len(proteins)

    def build_first_layer(self, proteins):
        # construct first layer (of proteins)
        layer = []
        for prot in proteins:
            prot_node = CTNode(prot)
            layer.append(prot_node)

        # pad with necessarily absent dummy variables so that the
        # number of variables is a power of 2; this is not the most
        # efficient method for this. because they are absent, they won't influence the
        # total sum, and thus Ds.
        for i in range(0, 2 ** self.log_length - len(proteins)):
            # this protein cannot be present, therefor set probability array to (0,1)
            layer.append(CTNode([np.array([1, 0])]))  # TODO change this order

        self.all_layers.append(layer)

    def build_remaining_layers(self):
        # construct layers of count nodes
        for L in range(self.log_length):
            # print('layers needed: ',int(len(self.allLayers[0])/(2**(L+1))))
            most_recent_layer = self.all_layers[-1]
            layer = []
            for i in range(int(len(self.all_layers[0]) / (2 ** (L + 1)))):
                left_parent = most_recent_layer[i * 2]
                right_parent = most_recent_layer[i * 2 + 1]
                count_node = CTNode.create_count_node(left_parent, right_parent)
                layer.append(count_node)

            # add connection to remaining nodes (when layer above is not a power of 2)
            self.all_layers.append(layer)

        # final node gets (Ds | N) multiplied into its likelihoodBelow
        final_node = self.all_layers[-1][0]
        # normalize for greater precision
        final_node.likelihood_below = array_utils.normalize(self.n_to_shared_likelihoods)
        self.last_node = final_node

    def propagate_backward(self):
        # propagate backward, setting likelihoodBelow.
        # the loop has upper bound at logLength+1
        # because of the layer of proteins
        for L in range(1, self.log_length + 1)[::-1]:
            layer = self.all_layers[L]

            for i in range(len(layer)):
                node = layer[i]
                left_parent = node.left_parent
                right_parent = node.right_parent

                left_parent.likelihood_below = node.message_up_left()
                right_parent.likelihood_below = node.message_up_right()

        self.protein_layer = self.all_layers[0]

    def posterior_for_variable(self, prot_idx):
        return self.protein_layer[prot_idx].posterior()

    def message_to_variable(self, prot_idx):
        return self.protein_layer[prot_idx].messages_up()

    def message_to_shared_likelihood(self):
        return self.last_node.joint_above[0][0: (self.n_proteins + 1)]
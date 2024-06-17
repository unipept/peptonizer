import numpy as np

from scipy.special import logsumexp


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

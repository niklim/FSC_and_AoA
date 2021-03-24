import multiprocessing as mp
from collections import defaultdict
from sklearn.metrics.pairwise import cosine_similarity


"""
This module maps words, encoded orthographically or phonologically, to their form-semantic consistency values: 
form-based neighbors are words which embed the target (as in Marelli & Amenta).
"""


embeddings = None
frequencies = None


def _mp_compute_te_fsc(args):

    return compute_te_fsc(*args)


def compute_te_fsc(target, neighbors):

    num, den = 0, 0
    target_vector = embeddings.get_vector(target)
    for neighbor in neighbors:
        den += frequencies[neighbor]
        try:
            neighbour_vector = embeddings.get_vector(neighbor)
            sim = abs(cosine_similarity(target_vector, neighbour_vector))
            num += sim * frequencies[neighbor]
        except ValueError:
            continue
    fsc = float(num / den)

    return target, fsc


def target_embedded_fsc(targets2neighbors, embedding_space, word2freq, threads=32):

    """
    :param targets2neighbors:   dict, mapping target words to lists of neighbors
    :param embedding_space:     SemanticSpace object containing the semantic embeddings
    :param word2freq:           dict, mapping target words to frequency values from SUBTLEX
    :param threads:             int, indicating how many cores to spread the process on
    :return:                    dict, target words mapped to form-semantic consistency values computed using
                                target-embedded neihbors, as in Marelli & Amenta (2015 and 2017)
    """

    global embeddings
    global frequencies
    embeddings = embedding_space
    frequencies = word2freq

    t2fsc = defaultdict(float)

    with mp.Pool(threads) as pool:
        outputs = pool.imap(
            _mp_compute_te_fsc, ((target, neighbors) for target, neighbors in targets2neighbors.items())
        )
        for output in outputs:
            t, fsc = output
            t2fsc[t] = fsc

    embeddings, filter_vocab = None, None

    return t2fsc

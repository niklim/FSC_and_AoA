import multiprocessing as mp
from collections import defaultdict
from sklearn.metrics.pairwise import cosine_similarity


"""
This module maps words, encoded orthographically or phonologically, to their form-semantic consistency values:
form-based neighbors are words with low Levenshtein distance to the target (as in Hendrix and Sun 2020.
"""


embeddings = None


def _mp_compute_ld_fsc(args):

    return compute_ld_fsc(*args)


def compute_ld_fsc(target, neighbors):

    num, den = 0, len(neighbors)
    target_vector = embeddings.get_vector(target)

    for n, lev_dis in neighbors:
        try:
            neighbour_vector = embeddings.get_vector(n)
            sim = abs(cosine_similarity(target_vector, neighbour_vector))
            num += sim / lev_dis
        except ValueError:
            continue
    fsc = float(num / den)

    return target, fsc


def levenshtein_fsc(targets2neighbors, embedding_space, threads=32):

    """
    :param targets2neighbors:   dict, mapping target words to lists of neighbors
    :param embedding_space:     SemanticSpace object containing the semantic embeddings
    :param threads:             int, the number of cores to spread the process over
    :return:                    dict, target words mapped to form-semantic consistency values computed using levenshtein
                                distance to find nearest neighbors based on form, as in Hendrix and Sun 2020
    """

    global embeddings
    embeddings = embedding_space

    t2fsc = defaultdict(int)

    with mp.Pool(threads) as pool:
        outputs = pool.imap(
            _mp_compute_ld_fsc, ((target, neighbors) for target, neighbors in targets2neighbors.items())
        )
        for output in outputs:
            t, fsc = output
            t2fsc[t] = fsc
    embeddings = None

    return t2fsc

import numpy as np
import multiprocessing as mp
from datetime import datetime
from collections import defaultdict

embeddings = None
filter_vocab = None
n_neighbors = None


def _mp_compute_snd(args):

    return compute_snd(args)


def compute_snd(w):

    distances_df = embeddings.all_distances([w])
    filtered_distances_df = distances_df.loc[:, distances_df.columns.isin(filter_vocab)]
    sorted_distances = filtered_distances_df.loc[w].sort_values(ascending=True)
    snd = np.mean(sorted_distances[1:n_neighbors + 1])  # index 0 is identical to target word, so skip it

    return w, snd


def neighborhood_density(sem_space, target_words, reference_vocab, n=20, threads=32):

    """
    :param sem_space:       a SemanticSpace object.
    :param target_words:    list, target words for which to compute snd.
    :param reference_vocab: list, reference vocabulary listing words to consider as valid neighbors.
    :param n:               int, number of neighbours to consider. Default to 20.
    :param threads:         int, the number of cores to use for parallel processing.
    :return:                dict, mapping words to their respective SND values. Higher values indicate sparser semantic
                            neighborhoods
    """

    global embeddings
    global filter_vocab
    global n_neighbors
    embeddings = sem_space
    filter_vocab = reference_vocab
    n_neighbors = n

    w2snd = defaultdict(float)

    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Started computing semantic neighborhood density..."))

    with mp.Pool(threads) as pool:
        outputs = pool.imap(_mp_compute_snd, ((word) for word in target_words))
        for output in outputs:
            w, snd = output
            w2snd[w] = snd
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done."))

    embeddings, filter_vocab, n_neighbors = None, None, None

    return w2snd

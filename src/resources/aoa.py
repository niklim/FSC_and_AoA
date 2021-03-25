import pandas as pd


def read(path):

    """
    :param path:        the path to the file AoA_ratings_from_all_sources.xlsx from the Age of Acquisition norms
                        collected by Kuperman et al 2012. The exact file can be downloaded here:
                        http://crr.ugent.be/papers/AoA_ratings_from_all_sources.zip
    :return:            the set of words for which age of acquisition norms were collected.
    """

    aoa_norms = pd.read_excel(path, usecols = ["Word", "Rating.Mean"])
    aoa_words = set(aoa_norms['Word'])

    return aoa_words, aoa_norms

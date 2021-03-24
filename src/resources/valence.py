import pandas as pd


def read(path):

    """
    :param path:            the path to the csv file containing the valence values from Brysbaert et al.
    :return:                a dict measures words to valence scores
    """

    df = pd.read_csv(path, header=0)
    df["word"] = df["Word"].str.lower()
    word2val = pd.Series(df['V.Mean.Sum'].values, index=df['word']).to_dict()

    return word2val

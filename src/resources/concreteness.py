import pandas as pd


def read(path):

    """
    :param path:            the path to the .txt file containing the concreteness norms from Brysbaert et al.
    :return:                a dict measures words to concreteness scores
    """

    df = pd.read_csv(path, sep='\t', header=0)
    df["word"] = df["Word"].str.lower()
    word2concr = pd.Series(df['Conc.M'].values, index=df['word']).to_dict()

    return word2concr

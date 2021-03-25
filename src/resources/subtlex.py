import pandas as pd


def read(path):

    """
    :param path:            the path to the csv file containing the SUBTLEX-US frequency values
    :return:                a dict measures words to frequency values
    """

    df = pd.read_csv(path, header=0)
    df["word"] = df["Word"].str.lower()
    word2freq = pd.Series(df['SUBTLWF'].values, index=df['word']).to_dict()

    return word2freq

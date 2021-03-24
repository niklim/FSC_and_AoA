import pandas as pd


def read(path):

    """
    :param path:            the path to the csv file containing the OLD20 values, with words in the first column and
                            OLD20 values in the second, separated using white spaces.
    :return:                dict, mapping words to OLD20 scores
    """

    df = pd.read_csv(path, sep=' ', header=None)
    word2old = pd.Series(df[1].values, index=df[0]).to_dict()

    return word2old

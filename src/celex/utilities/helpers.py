__author__ = 'GCassani'

"""Helper functions that hardcode some of the information necessary to build the dictionary from the
   Celex database and fine-tune it to the processing of CHILDES corpora"""


def inflection_dict():

    """
    :return inflections:    a dictionary measures Celex inflectional codes to their meanings
    """

    inflections = {'S': 'SINGULAR',
                   'P': 'PLURAL',
                   'i': 'INFINITIVE',
                   'p': 'PARTICIPLE',
                   'e': 'PRESENT',
                   'a': 'PAST',
                   '3': 'THIRDPERSON'}

    return inflections


########################################################################################################################


def vowels():

    """
    :return celex_vowels:   a set containing the symbols indicating vowels in the Celex phonetic transcription
    """

    celex_vowels = {'I', 'E', '{', 'V', 'Q', 'U', '@', 'i', '#', '$', 'u', '3',
                    '1', '2', '4', '5', '6', '7', '8', '9', 'c', 'q', 'O', '~'}

    return celex_vowels

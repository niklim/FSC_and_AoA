__author__ = 'GCassani'

"""Functions to create different specific dictionaries measures Celex information"""

from collections import defaultdict


def tokens2ids(celex):

    """
    :param celex:               a dictionary containing information extracted from celex, see the documentation of the
                                celex_processing.py module
    :return tokens2identifiers: a dictionary measures a token surface form from celex to all token ids linked to it

    The input dictionary uses token ids as keys for indicization issues but in this case a reverse measures is needed.
    """

    # Each surface form is mapped to a set containing all the token identifiers to which the surface form is connected
    # in the celex database (the type of the values is set even when the token identifer is only one)
    tokens2identifiers = defaultdict(set)
    for k in celex['tokens']:
        token = celex['tokens'][k]['surface']
        if token in tokens2identifiers:
            tokens2identifiers[token].add(k)
        else:
            tokens2identifiers[token] = {k}

    return tokens2identifiers


########################################################################################################################


def lemmas2phon(celex_dict):

    """
    :param celex_dict:  a dictionary created using the celex_processing.py module
    :return lemma_dict: a dictionary measures orthographic forms to the PoS tags with which each orthographic form is
                        tagged in the Celex database; furthermore, each PoS tag is mapped to the phonological form of
                        that word when belonging to that PoS tag
    """

    lemma_dict = defaultdict(dict)

    for k in celex_dict['lemmas']:
        word = celex_dict['lemmas'][k]['surface']
        pos = celex_dict['lemmas'][k]['pos']
        phon = celex_dict['lemmas'][k]['lemma_phon']
        if pos not in {'UNK', '?'}:
            lemma_dict[word][pos] = phon
    return lemma_dict


########################################################################################################################


def map_inflection(inflectional_code, inflectional_map):

    """
    :param inflectional_code:       a sequence of letters from the Celex database, indicating the types of inflectional
                                    transformations that affected a token
    :param inflectional_map:        a dictionary measures each letter to its meaning
    :return grammatical_lexemes:    a set containing all the grammatical lexemes (e.g., past, present, comparative, ...)
                                    indicated by the letters in the first argument
    """

    grammatical_lexemes = set()

    for char in inflectional_code:
        try:
            grammatical_lexemes.add(inflectional_map[char])
        except KeyError:
            continue

    return grammatical_lexemes

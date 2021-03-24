import re
import pandas as pd
from collections import defaultdict


def read_mono(path):

    """
    :param path:    str, indicating the path to the morpholex database
    :return:        a set of lower-cased strings indicating the mono-morphemic words found in the morpholex database
    """

    morpholex_df = pd.read_excel(path, sheet_name=1)
    w = set(morpholex_df['MorphoLexSegm'])
    targets = set()
    regex = re.compile('[^a-z]')
    for el in w:
        try:
            targets.add(regex.sub('', el.lower()))
        except AttributeError:
            pass

    return targets


def read_mono_inflected(path):

    """
    :param path:    str, indicating the path to the morpholex database
    :return:        a set of lower-cased strings indicating the mono-morphemic words with inflectional morphemes
                        found in the morpholex database
    """

    targets = set()
    inflections = ['s', 'ed', 'ing', 'en', "'s", 'er', 'est', 'es', 'ies', 'ings', 'ied']
    morpholex = pd.read_excel(path, sheet_name=None)
    for name, sheet in morpholex.items():
        if name == '0-1-0':
            token2base = pd.Series(sheet['MorphoLexSegm'].values, index=sheet['Word']).to_dict()
            for token, base in token2base.items():
                if type(token) == int or type(token) == float:
                    continue
                base = base.strip('{{()}}')
                token = token.lower()
                if token != base:
                    inflected = [''.join([base, affix]) for affix in inflections]
                    reduplicated_final = ''.join([token, token[-1]])
                    inflected.extend([''.join([reduplicated_final, affix]) for affix in inflections])
                    if base[-1] == 'f' or token[-2:] == 'fe':
                        f_v_alternation = re.sub('(f|fe)$', 'v', base)
                        inflected.extend([''.join([f_v_alternation, affix]) for affix in inflections])
                    if base[-1] == 'e':
                        no_e = re.sub('e$', '', base)
                        inflected.extend([''.join([no_e, affix]) for affix in inflections])
                    if base[-1] == 'y':
                        no_y = re.sub('y$', '', base)
                        inflected.extend([''.join([no_y, affix]) for affix in inflections])
                    if token == 'vertices':
                        inflected.append(token)
                    if base.endswith('c'):
                        plus_k = base + 'k'
                        inflected.extend([''.join([plus_k, affix]) for affix in inflections])
                    if base[-1] in {'b', 'd', 'f', 'g', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'z'}:
                        reduplicated = base + base[-1]
                        inflected.extend([''.join([reduplicated, affix]) for affix in inflections])
                    if base.endswith('eau'):
                        inflected.extend([''.join([base, affix]) for affix in ['s', 'x']])

                    o_to_ou = re.sub('o[bcdfgklmnpqrstvxz]+$', 'ou', base)
                    inflected.extend([''.join([o_to_ou, affix]) for affix in inflections])
                    if token in inflected:
                        targets.add(token)

    return targets


def read_poly(path):

    """
    :param path:    str, indicating the path to the morpholex database
    :return:        a set of lower-cased strings indicating the poly-morphemic words found in the morpholex database
    """

    targets = set()
    regex = re.compile('[^a-z]')
    morpholex = pd.read_excel(path, sheet_name=None)
    for name, sheet in morpholex.items():
        try:
            a, b, c = name.split('-')
            if a != '0' or c != '0':
                w = set(sheet['MorphoLexSegm'])
                for el in w:
                    try:
                        targets.add(regex.sub('', el.lower()))
                    except AttributeError:
                        pass
        except ValueError:
            pass

    return targets


def compute_morph_complexity(shared_words, mono_words, mono_inflected_words, poly_words):

    """
    :param shared_words:            set (or list), the shared words of all data sets
    :param mono_words:              list, monomorphemic words extracted from MorphoLEX
    :param mono_inflected_words:    list, monomorphemic words with inflections extracted from MorphoLEX
    :param poly_words:              list, polymorphmeic words extracted from MorphoLEX
    :return:                        dict, words mapped to morphological status (monomorphemic = 0; polymorphemic = 1)
    """

    word2morphstatus = defaultdict(int)
    for word in shared_words:
        if word in mono_words or word in mono_inflected_words:
            word2morphstatus[word] = 0
        elif word in poly_words:
            word2morphstatus[word] = 1
        else:
            continue

    return word2morphstatus
__author__ = 'GCassani'

"""Function to get word level phonological information from the Celex database"""

from src.celex.phonology.utilities import get_stressed_vowel


def read_epw(epw_path, celex_dict, vowel_set, reduced=True, compounds=False):

    """
    :param epw_path:    the path to the epw.cd file from the CELEX database
    :param celex_dict:  a dictionary built using make_celex_dict (see make_celex_dict for the details)
    :param vowel_set:   a set containing the vowel characters used in the CELEX database to encode phonological
                        representations: CELEX contains 4 different encodings, here the DISC one is used and thus the
                        vowel set must contain all relevant vowels as encoded in DISC format
    :param reduced:     a boolean. If True, the function tries to get the reduced phonological versions every
                        time it can, to be closer to the spoken language. However, since not every token in CELEX has a
                        corresponding reduced version, when this try fails, the function falls back to the canonical
                        version of the token. If False, the standard phonological form is retrieved directly.
    :param compounds:   a boolean. If true, all entries in Celex are considered; if False, entries which contain spaces
                        are discarded
    :return celex_dict: the updated version of the celex dictionary, containing token identifiers as keys, each mapped
                        to a dictionary containing information about the token surface form, the lemma identifier of the
                        corresponding lemma, the phonological form of the token and the stressed vowel.

    ------------------------------------------------------------------------------
    |  EPW COLUMN KEYS                                                           |
    ------------------------------------------------------------------------------
    line[0]  : TokenID
    line[1]  : token surface form
    line[3]  : LemmaID
    line[4]  : number of pronunciations available
    line[5]  : status of pronunciation, primary vs secondary
    line[6]  : DISC phonetic encoding, with stress marker
    line[8]  : syllabified CELEX transcription with brackets
    ------------------------------------------------------------------------------
    |  NOT ALWAYS AVAILABLE                                                      |
    ------------------------------------------------------------------------------
    line[9]  : second pronunciation status, primary vs secondary
    line[10] : second pronunciation DISC phonetic encoding, with stress marker
    line[12] : second pronunciation syllabified CELEX transcription with brackets
    """

    # resources the epw.cd file
    with open(epw_path, 'r+') as epw:
        for line in epw:
            records = line.strip().split('\\')

            # if compounds is set to False and a space is detected in the string, skip
            if ' ' in records[1] and not compounds:
                continue

            else:
                # store the token identifier as key in the sub-dictionary 'tokens'
                # each token identifier is mapped to a dictionary containing several keys to which values are mapped,
                # including the token surface form ('surface'), the lemma identifier ('lemmaID'), the phonological form
                # ('phon'), and the stressed vowel ('stressed_vowel'). First assign values to the token surface form and
                # the lemma identifier.
                celex_dict['tokens'][records[0]]['surface'] = records[1]
                celex_dict['tokens'][records[0]]['lemmaID'] = records[3]

                # if reduced forms are allowed, try to fetch it and fall back onto the regular one if there isn't;
                # otherwise, go directly for the regular form
                if reduced:
                    try:
                        celex_dict['tokens'][records[0]]['phon'] = records[10]
                    except IndexError:
                        celex_dict['tokens'][records[0]]['phon'] = records[6]
                else:
                    celex_dict['tokens'][records[0]]['phon'] = records[6]

                # get the stressed vowel from the token phonological representation
                stressed_vowel = get_stressed_vowel(celex_dict['tokens'][records[0]]['phon'], vowel_set)
                celex_dict['tokens'][records[0]]['stressed_vowel'] = stressed_vowel

    return celex_dict

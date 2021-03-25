__author__ = 'GCassani'

"""Helper functions to get morphological information from Celex entries"""

import re


def get_constituent_morphemes(morphological_analyses):

    """
    :param morphological_analyses:  a string extracted from the eml.cd file from the CELEX database, which contains the
                                    morphological segmentation of a lemma, with brackets signaling the correct
                                    composition and a Part-of-Speech tag in brackets at the end of each morpheme. The
                                    outermost tag refers to the whole lemma.
    :return morphemes_clean:        a list (order is preserved) of the morphemes that make up the lemma.

    This function gets a complete but flat morphological analyses from the structured one that is retrieved from the
    CELEX database: the goal is to have all the consituent morphemes, but the composition scheme is not important.
    """

    morphemes_clean = []

    # get rid of parantheses by getting all the strings of lowercase letters (PoS tags are capitalized, and we don't
    # want them here)
    morphemes = re.findall(r"\(([a-z]+)\)", morphological_analyses)

    # given that the morphological analyses includes the full composition scheme, morphemes may be repeated: only add a
    # morpheme to the leave_k_out_mapping structure if the same morpheme isn't already included therein.
    for morpheme in morphemes:
        if morpheme not in morphemes_clean:
            morphemes_clean.append(morpheme)

    return morphemes_clean


########################################################################################################################


def get_part_of_speech(morphological_analyses):

    """
    :param morphological_analyses:  a string extracted from the eml.cd file from the CELEX database, which contains the
                                    morphological segmentation of a lemma, with brackets signaling the correct
                                    composition and a Part-of-Speech tag in brackets at the end of each morpheme. The
                                    outermost tag refers to the whole lemma.
    :return pos_tag:                a capital letter marking the lexical category of a lemma

    This function simply returns the PoS tag given the morphological analyses of a lemma extracted from the eml.cd file
    from the CELEX database.
    """

    # if no PoS tag is provided, return the label 'UNK'
    try:
        pos_tag = morphological_analyses[-2]
    except IndexError:
        pos_tag = 'UNK'

    return pos_tag

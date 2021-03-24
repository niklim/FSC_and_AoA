__author__ = 'GCassani'

"""Functions to get lemma level morphological information  from the Celex dictionary"""

from src.celex.morphology.utilities import get_constituent_morphemes, get_part_of_speech


def read_eml(eml_path, celex_dict):

    """
    :param eml_path:    the path to the eml.cd file from the CELEX database
    :param celex_dict:  a dictionary built using make_celex_dict and already updated using read_epw(), read_epl(),
                        and read_emw()
    :return celex_dict: the input dictionary, updated with information about the lemma that corresponds to every token
                        and its Part-of-Speech


    ------------------------------------------------------------------------------
    |  EMW COLUMN KEYS                                                           |
    ------------------------------------------------------------------------------
    line[0]  : LemmaID
    line[1]  : lemma surface form
    line[3]  : morphological status
                  'C' -> morphologically complex
                  'M' -> monomorphemic
                  'Z' -> zero derivation
                  'F' -> contracted form
                  'I' -> irrelevant morphology
                  'O' -> obscure morphology
                  'R' -> may include a root
                  'U' -> undetermined
    line[21]  : structured segmentation with brackets and PoS tag of the first,
                  default parsing of the lemma:
                  e.g. revolution = ((revolt),(ution)[N|V.])[N]
                  where the final [N] tells revolution is a noun and the
                  bracketing tells that it consists of two morphomes, the
                  verb revolt and the affix ution, which is specified to be
                  a suffix to the verb.
    line[40, 59, +19...] contain structured segmentations with brackets and Pos
                         tag of alternative parsings of the lemma
    """

    # collect all lemma identifiers that are already present in the celex dictionary
    lemma_identifiers = set()
    for k in celex_dict['tokens']:
        lemma_identifiers.add(celex_dict['tokens'][k]['lemmaID'])

    # resources the eml.cd file
    with open(eml_path, 'r+') as eml:
        for line in eml:
            records = line.strip().split('\\')

            # check if the lemma being considered exists in the celex dictionary
            # if it does, get its morphological analysis, and use it to derive the lemma PoS and its constituent
            # morphemes. Store these pieces of information in the celex dictionary
            if records[0] in lemma_identifiers:

                morphological_analyses = records[21]

                pos_tag = get_part_of_speech(morphological_analyses)
                morphemes = get_constituent_morphemes(morphological_analyses)

                celex_dict['lemmas'][records[0]]['surface'] = records[1]
                celex_dict['lemmas'][records[0]]['morph_analysis'] = morphemes
                celex_dict['lemmas'][records[0]]['pos'] = pos_tag

            # if the lemma being considered doesn't exist in the celex dictionary, keep running
            else:
                continue

    return celex_dict

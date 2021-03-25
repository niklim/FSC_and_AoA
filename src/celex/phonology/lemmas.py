__author__ = 'GCassani'

"""Function to extract lemma level phonological information from Celex"""


def read_epl(epl_path, celex_dict):

    """
    :param epl_path:    the path to the emw.cd file from the CELEX database
    :param celex_dict:  a dictionary built using make_celex_dict and already updated using read_epw()
    :return celex_dict: the input dictionary, updated with information about the phonological form of each lemma
    """

    lemma_identifiers = set()
    for k in celex_dict['tokens']:
        lemma_identifiers.add(celex_dict['tokens'][k]['lemmaID'])

    with open(epl_path, 'r+') as epl:
        for line in epl:
            records = line.strip().split('\\')

            if records[0] in lemma_identifiers:
                celex_dict['lemmas'][records[0]]['lemma_phon'] = records[5]

    return celex_dict

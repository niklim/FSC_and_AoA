__author__ = 'GCassani'

"""Function to get word level morphological information from the Celex database"""


def read_emw(emw_path, celex_dict, compounds=False):

    """
    :param emw_path:    the path to the emw.cd file from the CELEX database
    :param celex_dict:  a dictionary built using make_celex_dict and already updated using read_epw() and read_epl()
    :param compounds:   a boolean. If true, all entries in Celex are considered; if False, entries which contain spaces
                        are discarded
    :return celex_dict: the input dictionary, updated with information about the inflectional morphology of each token

    ------------------------------------------------------------------------------
    |  EMW COLUMN KEYS                                                           |
    ------------------------------------------------------------------------------
    line[0]  : TokenID
    line[1]  : token surface form
    line[3]  : LemmaID
    line[4]  : flectional type
                  'S' -> singular
                  'P' -> plural
                  'b' -> positive
                  'c' -> comparative
                  's' -> superlative
                  'i' -> infinitive
                  'p' -> participle
                  'e' -> present
                  'a' -> past
                  '1' -> first person
                  '2' -> second person
                  '3' -> third person
                  'r' -> rare
                  'X' -> head (no N, V, Adj, Adv)
    line[5]  : flectional variant, which letters are removed or added to form the
                  inflected form [-(.) for deleted letters, +(.) for added letters,
                  @ for the stem]
    """

    # collect all token identifiers already stored in the celex dictionary
    token_ids = set()
    for k in celex_dict['tokens']:
        token_ids.add(k)

    # resources the emw.cd file
    with open(emw_path, 'r+') as emw:
        for line in emw:
            records = line.strip().split('\\')

            # if compounds is set to False and a space is detected in the string, skip
            if ' ' in records[1] and not compounds:
                continue
            else:
                # check that the token being considered already exists in the celex dictionary and then that its surface
                # form and lemma identifier are the same as those of the corresponding entry in the celex dictionary. If
                # all these conditions hold, get the inflection codes associated with the token and store them in the
                # celex dictionary
                if records[0] in token_ids:
                    if records[1] == celex_dict['tokens'][records[0]]['surface'] \
                            and records[3] == celex_dict['tokens'][records[0]]['lemmaID']:
                        celex_dict['tokens'][records[0]]['inflection'] = records[4]
                    else:
                        print("TokenID: " + records[0] +
                              " is in celex dict but the surface form or the lemmaID don't match.")
                else:
                    print("TokenID: " + records[0] + " couldn't be located in celex dict.")

    return celex_dict

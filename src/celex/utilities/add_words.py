__author__ = 'GCassani'

"""This functions hard-codes extra words that are not in Celex , to increase the coverage of Celex"""


def hardcode_words(celex_dict, token_surface, token_id, lemma_id, token_phonetic, vowel, inflection, lemma_surface,
                   lemma_phonetic, morph, pos):

    """
    :param celex_dict:      a dictionary created using the function initialize_celex_dict from this module. The
                            dictionary can be empty or already contain items
    :param token_surface:   a string containing the orthographic form of the token you want to hard-code
    :param token_id:        the numeric ID (encoded as string) of the token to be added. Make sure you choose an id
                            that doesn't already exist in CELEX
    :param lemma_id:        the numeric ID (encoded as string) of the lemma to which the token to be hard-coded
                            corresponds. If the lemma already exists in CELEX, make sure you use the correct ID
    :param token_phonetic:  a string containing the phonetic form of the token to be hard-coded - make sure you use the
                            charset used in CELEX
    :param vowel:           a string containing the vowel bearing the stress in the token to be hard-coded
    :param inflection:      a string indicating the inflectional code(s) of the word to be hard-coded
    :param lemma_surface:   the orthographic form of the corresponding lemma
    :param lemma_phonetic:  the phonetic form of the corresponding lemma
    :param morph:           a list of strings indicating the morphological analysis of the lemma
    :param pos:             the Part-of-Speech tag of the lemma
    :return celex_dict:     the dictionary updated with the word to be hard-coded
    """

    celex_dict['tokens'][token_id]['surface'] = token_surface
    celex_dict['tokens'][token_id]['lemmaID'] = lemma_id
    celex_dict['tokens'][token_id]['phon'] = token_phonetic
    celex_dict['tokens'][token_id]['stressed_vowel'] = vowel
    celex_dict['tokens'][token_id]['inflection'] = inflection
    celex_dict['lemmas'][lemma_id]['surface'] = lemma_surface
    celex_dict['lemmas'][lemma_id]['morph_analysis'] = morph
    celex_dict['lemmas'][lemma_id]['lemma_phon'] = lemma_phonetic
    celex_dict['lemmas'][lemma_id]['pos'] = pos

    return celex_dict


########################################################################################################################


def add_childes_words(celex_dict):

    """
    :param celex_dict:  the Celex dictionary (can be empty)
    :return celex_dict: the input dictionary, to which several hardcoded words were added to increase coverage for
                        CHILDES transcripts
    """

    celex_dict = hardcode_words(celex_dict, "will'nt", '1000000', '51924', "'wIl-Ht", 'I', 'X', "won't",
                                "'w5nt", ['will', 'not'], 'V')
    celex_dict = hardcode_words(celex_dict, "mhm", '1000001', '1000001', "em", '_', 'X', "non_ling",
                                '_', ['non_ling'], 'C')
    celex_dict = hardcode_words(celex_dict, "lego", '1000002', '1000002', "'l3-go", 'e', 'S', "lego",
                                "'l3-go", ['lego'], 'N')
    celex_dict = hardcode_words(celex_dict, "colour", '1000003', '8425', "'kV-l@R", 'V', 'S', "colour",
                                "'kV-l@R", ['colour'], 'N')
    celex_dict = hardcode_words(celex_dict, "color", '1000004', '8425', "'kV-l@R", 'V', 'S', "colour",
                                "'kV-l@R", ['colour'], 'N')
    celex_dict = hardcode_words(celex_dict, "horsie", '1000005', '21619', "h$s", '$', 'S', "horse",
                                "'h$s", ['horse'], 'N')
    celex_dict = hardcode_words(celex_dict, "ssh", '1000006', '1000001', "C", '_', 'X', "non-ling",
                                '_', ['non-ling'], 'C')
    celex_dict = hardcode_words(celex_dict, "byebye", '1000007', '5863', "b2-'b2", '2', 'X', "bye-bye",
                                "b2-'b2", ['bye'], 'C')
    celex_dict = hardcode_words(celex_dict, "shooter", '1000008', '1000003', "'Su-t@R", 'u', 'S', "shooter",
                                "'Su-t@R", ['shoot', 'er'], 'N')
    celex_dict = hardcode_words(celex_dict, "mousie", '1000009', '29286', "'m2-sI", '2', 'S', "mouse",
                                "'m6s", ['mouse', 'y'], 'N')
    celex_dict = hardcode_words(celex_dict, "mummie", '1000010', '29460', "'mV-mI", 'V', 'S', "mummy",
                                "'mV-mI", ['mum', 'y'], 'N')
    celex_dict = hardcode_words(celex_dict, "favorite", '1000011', '16271', "'f1-vrIt", '1', 'b', "favourite",
                                "'f1-v@-rIt", ['favour', 'ite'], 'A')
    celex_dict = hardcode_words(celex_dict, "upside", '1000012', '1000004', "Vp-'s2d", '2', 'b', "upside",
                                "Vp-'s2d", ['upside'], 'B')
    celex_dict = hardcode_words(celex_dict, "carwash", '1000013', '1000005', "'k#R-wQS", '#', 'S', "carwash",
                                "'k#R-wQS", ['car', 'wash'], 'N')
    celex_dict = hardcode_words(celex_dict, "anymore", '1000014', '1000006', "E-nI-'m$R", '$', 'X', "anymore",
                                "E-nI-'m$R", ['any', 'more'], 'B')
    celex_dict = hardcode_words(celex_dict, "whee", '1000015', '1000001', "wee", '_', 'X', "non-ling",
                                '_', ['non-ling'], 'C')
    celex_dict = hardcode_words(celex_dict, "carpark", '1000016', '1000007', "k#R-'p#k", '#', 'S', "carpark",
                                "k#R-'p#k", ['car', 'park'], 'N')
    celex_dict = hardcode_words(celex_dict, "lawnmower", '1000017', '1000008', "'l$n-m5-er", '$', 'S', "lawnmower",
                                "'l$n-m5-er", ['lawn', 'mow', 'er'], 'N')
    celex_dict = hardcode_words(celex_dict, "whoo", '1000018', '1000001', "hU", '_', 'X', "non-ling",
                                '_', ['non-ling'], 'C')
    celex_dict = hardcode_words(celex_dict, "doggie", '1000019', '13205', "'dQ-gI", 'Q', 'S', "doggy",
                                "'dQ-gI", ['dog', 'y'], 'N')
    celex_dict = hardcode_words(celex_dict, "hotdog", '1000020', '21690', "hQt-'dQg", 'Q', 'S', "hot dog",
                                "hQt-'dQg", ['hot', 'dog'], 'N')
    celex_dict = hardcode_words(celex_dict, "christmas", '1000021', '7479', "'krIs-m@s", 'I', 'X', "christmas",
                                "'krIs-m@s", ['christ', 'mas'], 'N')
    celex_dict = hardcode_words(celex_dict, "traveling", '1000022', '48164', "'tr{v-lIN", '{', 'pe', "travel",
                                "'tr{-vP", ['travel'], 'V')
    celex_dict = hardcode_words(celex_dict, "snowplow", '1000023', '43060', "'snO-pl8", 'O', 'X', "snowplough",
                                "'sn5-pl6", ['snow', 'plough'], 'N')
    celex_dict = hardcode_words(celex_dict, "colors", '1000024', '8425', "'kV-l@Rs", 'V', 'P', "colour",
                                "'kV-l@R", ['colour'], 'N')

    return celex_dict

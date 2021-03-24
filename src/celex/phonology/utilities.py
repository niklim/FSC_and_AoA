__author__ = 'GCassani'

"""Helper function to extract phonological information from Celex entries"""


def get_stressed_vowel(phonetic_form, vowel_set):

    """
    :param phonetic_form:   a string containing the phonetic transcription of a token as extracted from the epw.cd file
                            from the CELEX database, using the DISC format and including stress markers. In this
                            encoding, the token is syllabified with syllables being divided by a dash ('-')
    :param vowel_set:          a set of vowel symbols used in the DISC encoding
    :return stressed_vowel: the vowel in the syllable bearing the stress marker (an apostrophe, " ' ")

    This function extracts the stressed vowel from the phonological representation of a token, as extracted from the
    CELEX database.
    """

    # get all the syllables that make up the full representation
    syllable = phonetic_form.split('-')

    # scan each syllable, check whether it contains the stress marker ("'"): if it does, get the vowel from the syllable
    # and return it. Each syllable only have one vowel, so no chance of picking the wrong one; moreover, no word has
    # more than one stressed syllable (primary stress, which is what we are interested in here). If a token doesn't
    # contain stressed vowels, return a dot ('.')
    for idx, syl in enumerate(syllable):
        if "'" in syl:
            for phoneme in syl:
                if phoneme in vowel_set:
                    stressed_vowel = phoneme
                    return stressed_vowel

    return '.'

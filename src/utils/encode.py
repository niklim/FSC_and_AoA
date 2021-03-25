from src.celex.utilities.dictionaries import tokens2ids


def get_celex_coverage(words, celex_dict):

    """
    :param words:       a set of target words to be encoded phonologically
    :param celex_dict:  the dictionary derived from the CELEX database
    :return:    - a set containing the words among the targets for which a unique phonological transcription was found
                - a set containing the words among the targets for which more than one phonological transcription was
                    found
                - a set containing the words among the targets for which no phonological transcription was found
    """

    tokens2identifiers = tokens2ids(celex_dict)

    no_phon = set()
    ambiguous = set()
    spoken_words_phonology = set()

    for word in words:
        token_ids = tokens2identifiers[word]
        possible_phonological_transcriptions = set()
        if token_ids:
            for token_id in token_ids:
                possible_phonological_transcriptions.add(celex_dict['tokens'][token_id]['phon'])
            if len(possible_phonological_transcriptions) > 1:
                ambiguous.add((word, tuple(possible_phonological_transcriptions)))
            else:
                for transcription in possible_phonological_transcriptions:
                    spoken_words_phonology.add((word, transcription.replace("-", "")))
        else:
            no_phon.add(word)

    return spoken_words_phonology, ambiguous, no_phon
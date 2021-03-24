import pandas as pd


def write_df(
        targets, out_path, d_aoa, d_OSC_te, d_OSC_ld, d_PSC_te, d_PSC_ld, d_conc, d_valence,
        d_freq_subtlex, d_word2old, d_target_phon, d_morph, d_snd
):

    """
    :param targets:         iterable, containing target words
    :param out_path:        str, indicating the path where to save the combined df
    :param d_aoa:           dict, maps target words to aoa values
    :param d_OSC_te:        dict, maps target words to OSC values computed using target embedded neighbors
    :param d_OSC_ld:        dict, maps target words to OSC values computed using levenshtein distance neighbors
    :param d_PSC_te:        dict, maps target words to PSC values computed using target embedded neighbors
    :param d_PSC_ld:        dict, maps target words to PSC values computed using levenshtein distance neighbors
    :param d_conc:          dict, maps target words to concreteness scores
    :param d_valence:       dict, maps target words to valence norms
    :param d_freq_subtlex:  dict, maps target words to frequency counts from SUBTLEX
    :param d_word2old:      dict, maps target words to OLD20 values
    :param d_target_phon:   dict, maps target words to corresponding phonological forms
    :param d_morph:         dict, maps target words to morphological status
    :param d_snd:           dict, maps target words to semantic neighborhood density values
    """
    values = []
    for word in targets:
        values.append(
            [word, len(word), len(d_target_phon[word]), d_aoa[word], d_OSC_te[word], d_OSC_ld[word], d_PSC_te[word],
             d_PSC_ld[word], d_conc[word], d_valence[word], d_freq_subtlex[word], d_word2old[word], d_snd[word],
             d_morph[word]]
        )

    final_df = pd.DataFrame(
        data=values,
        columns=[
            "word", "length", "n_phon", "aoa", "OSC_te", "OSC_ld", "PSC_te",
            "PSC_ld", "concr", "val", "freq", "old20", "snd", "morph"
        ]
    )

    final_df.to_csv(out_path, index=False)

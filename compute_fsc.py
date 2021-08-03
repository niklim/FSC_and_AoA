import os
import copy
import json
import random
import pickle
import itertools
import numpy as np
import pandas as pd
from datetime import datetime
from src.utils.io_ops import write_df
from semspaces.space import SemanticSpace
from src.utils.set_ops import shared_words
from src.utils.encode import get_celex_coverage
from src.measures.fsc_ld import levenshtein_fsc
from src.measures.fsc_te import target_embedded_fsc
from src.measures.semantic import neighborhood_density
from src.resources import aoa, concreteness, valence, subtlex, old20, morpholex
from src.utils.neighbors import get_target_embedded_neighbors, get_levenshtein_neighbours


random_baseline = False
sample_vocab = True

data_dir = os.path.join(os.getcwd(), 'data/')
out_dir = os.path.join(os.getcwd(), 'output/')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# LOAD DATA
print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Started loading the data..."))

embedding_space = SemanticSpace.from_csv(
    os.path.join(data_dir, '/home/gcassani/Resources/Embeddings/embedding_space.cbow.ukwac.subtlex.300dims.w5.w2v'),
    prenorm=True
)
w2v_words = embedding_space.included_words()
celex = json.load(open(os.path.join(data_dir, 'celex_dict.json')))

aoa_words, aoa_norms = aoa.read(os.path.join(data_dir, "AoA.xlsx"))
w2aoa = pd.Series(aoa_norms["Rating.Mean"].values, index=aoa_norms["Word"]).to_dict()

w2concr = concreteness.read(os.path.join(data_dir, "concreteness.txt"))
w2val = valence.read(os.path.join(data_dir, "valence.csv"))
w2freq = subtlex.read(os.path.join(data_dir, "subtlex.csv"))
w2old = old20.read(os.path.join(data_dir, "word2old.csv"))

mono = list(morpholex.read_mono(os.path.join(data_dir, "MorphoLEX_en.xlsx")))
poly = list(morpholex.read_poly(os.path.join(data_dir, "MorphoLEX_en.xlsx")))
mono_inflected = list(morpholex.read_mono_inflected(os.path.join(data_dir, "MorphoLEX_en.xlsx")))
print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done"))


# CONCATENATE THE DIFFERENT SETS OF POTENTIAL TARGET WORDS FOR LATER FILTERING STEPS
print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Started finding the target words..."))
filter_words = list(
    itertools.chain(
        aoa_words, set(w2concr.keys()), set(w2val.keys()), mono, poly, mono_inflected,
    )
)
print('The reference vocabulary for measuring SND consists of {} words.'.format(len(filter_words)))

# FIND WORDS SHARED ACROSS RESOURCES, SUCH THAT WE CAN ESTIMATE ALL NECESSARY PREDICTORS: FREQUENCY, CONCRETENESS,
# VALENCE, SEMANTIC NEIGHBOURHOOD DENSITY, AND OLD20. SINCE ICONICITY NORMS ARE FAR SMALLER IN SIZE, A SEPARATE ANALYSIS
# IS RUN WITH THOSE, WITHOUT FILTERING WORDS FOR THE MAIN EXPERIMENT
shared = shared_words(
    aoa_words, set(w2concr.keys()), set(w2val.keys()), list(w2v_words), set(w2freq.keys()), set(w2old.keys())
)
print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done"))

# ASSIGN BOOLEAN VALUE TO EACH WORD FOR WHICH ALL VARIABLES ARE AVAILABLE DEPENDING ON THEIR MORPHOLOGICAL COMPLEXITY
print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Mapping targets to psycholinguistic properties..."))
w2morph = morpholex.compute_morph_complexity(shared, mono, mono_inflected, poly)

# RESTRICT TO WORDS FOR WHICH A UNIQUE, UNAMBIGUOUS PHONOLOGICAL TRANSCRIPTION IS AVAILABLE IN CELEX, THEN MAP
# (ORTHO, PHONO) TUPLES TO AOA VALUES
t2phon = {k: v for k, v in get_celex_coverage(w2morph.keys(), celex)[0]}
t2aoa = {k: v for k, v in w2aoa.items() if k in t2phon.keys()}
print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done."))

# COMPUTE SEMANTIC NEIGHBORHOOD DENSITY
snd_path = os.path.join(data_dir, "target2snd.json")
try:
    t2snd = json.load(open(snd_path, "rb"))
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: File {} found and loaded.".format(snd_path)))
except FileNotFoundError:
    t2snd = neighborhood_density(embedding_space, t2phon.keys(), filter_words)
    json.dump(t2snd, open(snd_path, 'w'))


# FIND NEIGHBORS (TARGET-EMBEDDING AND LEVENSHTEIN) FOR ORTHOGRAPHIC AND PHONOLOGICAL FORMS
reference_vocab = [x for x in set(w2freq.keys()) if str(x) != 'nan']

print('The reference vocabulary for retrieving form-based neighbors consists of {} words.'.format(len(reference_vocab)))

ortho2neighbors_te = get_target_embedded_neighbors(set(t2phon.keys()), reference_vocab)
print(
    datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done retrieving target-embedded neighbors for orthographic neighbors.")
)
ortho2neighbors_ld = get_levenshtein_neighbours(set(t2phon.keys()), reference_vocab)
print(
    datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done retrieving Levenshtein distance neighbors for orthographic forms.")
)

phono2neighbors_te = get_target_embedded_neighbors(t2phon, reference_vocab, celex=celex)
print(
    datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done retrieving target-embedded neighbors for phonological forms.")
)
phono2neighbors_ld = get_levenshtein_neighbours(t2phon, reference_vocab, celex=celex)
print(
    datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done retrieving Levenshtein distance neighbors for phonological forms.")
)

# TRY TO FETCH FSC VALUES (ORTHO AND PHONO, COMPUTED USING TARGET-EMBEDDING NEIGHBORS OR LEVENSHTEIN DISTANCE
# NEIGHBORS) FROM FILE. IF THE FILE DOESN'T EXIST, COMPUTE VALUES
fsc_dir = os.path.join(out_dir, 'FSCmeasures')
if not os.path.exists(fsc_dir):
    os.makedirs(fsc_dir)

osc_te_path = os.path.join(fsc_dir, "OSC_te.pkl")
try:
    t2OSC_te = pickle.load(open(osc_te_path, "rb"))
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: File {} found and loaded.".format(snd_path)))
except FileNotFoundError:
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Started computing target-embedded OSC..."))
    t2OSC_te = target_embedded_fsc(ortho2neighbors_te, embedding_space, w2freq)
    pickle.dump(t2OSC_te, open(osc_te_path, "wb"))
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done."))
    print()

osc_ld_path = os.path.join(fsc_dir, "OSC_ld.pkl")
try:
    t2OSC_ld = pickle.load(open(osc_ld_path, "rb"))
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: File {} found and loaded.".format(snd_path)))
except FileNotFoundError:
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Started computing Levenshtein OSC..."))
    t2OSC_ld = levenshtein_fsc(ortho2neighbors_ld, embedding_space)
    pickle.dump(t2OSC_ld, open(osc_ld_path, "wb"))
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done."))
    print()

psc_te_path = os.path.join(fsc_dir, "PSC_te.pkl")
try:
    t2PSC_te = pickle.load(open(psc_te_path, "rb"))
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: File {} found and loaded.".format(snd_path)))
except FileNotFoundError:
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Started computing target-embedded PSC..."))
    t2PSC_te = target_embedded_fsc(phono2neighbors_te, embedding_space, w2freq)
    pickle.dump(t2PSC_te, open(psc_te_path, "wb"))
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done."))
    print()

psc_ld_path = os.path.join(fsc_dir, "PSC_ld.pkl")
try:
    t2PSC_ld = pickle.load(open(psc_ld_path, "rb"))
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: File {} found and loaded.".format(snd_path)))
except FileNotFoundError:
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Started computing Levenshtein PSC..."))
    t2PSC_ld = levenshtein_fsc(phono2neighbors_ld, embedding_space)
    pickle.dump(t2PSC_ld, open(psc_ld_path, "wb"))
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Done."))
    print()

# WRITE MEASURES TO FILE FOR SUBSEQUENT STATISTICAL ANALYSIS
write_df(
    t2phon.keys(), os.path.join(fsc_dir, "fsc_measures.csv"), w2aoa, t2OSC_te, t2OSC_ld, t2PSC_te,
    t2PSC_ld, w2concr, w2val, w2freq, w2old, t2phon, w2morph, t2snd
)

if random_baseline:

    n_subsamples = 1000
    seeds = random.sample(range(0, 100000000), n_subsamples)
    print(datetime.now().strftime(
        "%d/%m/%Y %H:%M:%S: Started computing FSC from {} random permutations of the embeddings...".format(n_subsamples)
    ))

    # COMPUTE FSC MEASURES FROM RANDOM PERMUTATIONS OF THE WORD EMBEDDINGS, REPEAT 1000 TIMES AND SAVE MEASURES TO FILE
    random_embeddings = copy.deepcopy(embedding_space)
    fsc_dir_rnd = os.path.join(out_dir, 'FSCrandom')
    if not os.path.exists(fsc_dir_rnd):
        os.makedirs(fsc_dir_rnd)

    for i, seed in enumerate(seeds):
        np.random.seed(seed)
        print(datetime.now().strftime("%d/%m/%Y %H:%M:%S: Started permutation {} of {}...".format(i + 1, n_subsamples)))
        random_embeddings.vectors = np.random.permutation(random_embeddings.vectors)

        t2OSC_te_rnd = target_embedded_fsc(ortho2neighbors_te, random_embeddings, w2freq)
        t2OSC_ld_rnd = levenshtein_fsc(ortho2neighbors_ld, random_embeddings)
        t2PSC_te_rnd = target_embedded_fsc(phono2neighbors_te, random_embeddings, w2freq)
        t2PSC_ld_rnd = levenshtein_fsc(phono2neighbors_ld, random_embeddings)

        write_df(
            t2phon.keys(), os.path.join(fsc_dir_rnd, "df{}.csv".format(i+1)), w2aoa, t2OSC_te_rnd, t2OSC_ld_rnd,
            t2PSC_te_rnd, t2PSC_ld_rnd, w2concr, w2val, w2freq, w2old, t2phon, w2morph, t2snd
        )


if sample_vocab:

    n_subsamples = 500
    sampling_rates = [50, 75]
    seeds = random.sample(range(0, 100000000), n_subsamples)

    fsc_dir_subset = os.path.join(out_dir, 'FSCsubset')
    if not os.path.exists(fsc_dir_subset):
        os.makedirs(fsc_dir_subset)

    for rate in sampling_rates:
        print(datetime.now().strftime(
            "%d/%m/%Y %H:%M:%S: Started computing FSC from {} samples of {}% of the reference vocabulary...".format(
                n_subsamples, rate
            )
        ))

        fsc_subdir_subset = os.path.join(fsc_dir_subset, 'rate{}'.format(rate))
        if not os.path.exists(fsc_subdir_subset):
            os.makedirs(fsc_subdir_subset)

        # determine the size of each random sample of the reference vocabulary given the sampling rate
        target_vocab_size = round(len(reference_vocab)*rate/100)

        for i, seed in enumerate(seeds):
            np.random.seed(seed)
            print(datetime.now().strftime(
                "%d/%m/%Y %H:%M:%S: Sample {} of {}...".format(i + 1, n_subsamples))
            )
            reference_vocab_subsample = random.sample(reference_vocab, target_vocab_size)

            ortho2neighbors_te_subsample = get_target_embedded_neighbors(set(t2phon.keys()), reference_vocab_subsample)
            print(
                datetime.now().strftime(
                    "%d/%m/%Y %H:%M:%S: Done retrieving target-embedded neighbors for orthographic neighbors.")
            )
            ortho2neighbors_ld_subsample = get_levenshtein_neighbours(set(t2phon.keys()), reference_vocab_subsample)
            print(
                datetime.now().strftime(
                    "%d/%m/%Y %H:%M:%S: Done retrieving Levenshtein distance neighbors for orthographic forms.")
            )

            phono2neighbors_te_subsample = get_target_embedded_neighbors(t2phon, reference_vocab_subsample, celex=celex)
            print(
                datetime.now().strftime(
                    "%d/%m/%Y %H:%M:%S: Done retrieving target-embedded neighbors for phonological forms.")
            )
            phono2neighbors_ld_subsample = get_levenshtein_neighbours(t2phon, reference_vocab_subsample, celex=celex)
            print(
                datetime.now().strftime(
                    "%d/%m/%Y %H:%M:%S: Done retrieving Levenshtein distance neighbors for phonological forms.")
            )

            t2OSC_te_subset = target_embedded_fsc(ortho2neighbors_te_subsample, embedding_space, w2freq)
            t2OSC_ld_subset = levenshtein_fsc(ortho2neighbors_ld_subsample, embedding_space)
            t2PSC_te_subset = target_embedded_fsc(phono2neighbors_te_subsample, embedding_space, w2freq)
            t2PSC_ld_subset = levenshtein_fsc(phono2neighbors_ld_subsample, embedding_space)

            write_df(
                t2phon.keys(), os.path.join(fsc_subdir_subset, "df{}.csv".format(i + 1)), w2aoa, t2OSC_te_subset,
                t2OSC_ld_subset, t2PSC_te_subset, t2PSC_ld_subset, w2concr, w2val, w2freq, w2old, t2phon, w2morph, t2snd
            )


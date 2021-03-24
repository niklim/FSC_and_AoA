__author__ = 'GCassani'

"""Functions to initialize, create, print, and load the Celex dictionary"""

import os
import json
from time import strftime
from collections import defaultdict
from src.celex.phonology.words import read_epw
from src.celex.phonology.lemmas import read_epl
from src.celex.morphology.words import read_emw
from src.celex.morphology.lemmas import read_eml
from src.celex.utilities.helpers import vowels
from src.celex.utilities.add_words import add_childes_words


def get_name(compounds=False, reduced=False):

    dict_name = 'celex_dict'
    if compounds:
        dict_name += '_compounds'
    if reduced:
        dict_name += '_reduced'

    return '.'.join([dict_name, 'json'])


########################################################################################################################


def initialize_celex_dict():

    """
    :return celex_dict: this function makes a dictionary consisting of two dictionaries, 'tokens' and 'lemmas', which
                        in turn consists of dictionaries of dictionaries.
    """

    celex_dict = {'tokens': defaultdict(dict),
                  'lemmas': defaultdict(dict)}
    return celex_dict


########################################################################################################################


def print_celex_dict(celex_dict, outfile):

    """
    :param celex_dict:  the dictionary created using make_celex_dictionary(). It can or cannot have been updated using
                        the functions read_epw, read_emw, and read_eml: this function always tries to print information
                        from the dictionary but whenever it cannot find it, it substitutes the missing information with
                        a dash ('-')
    :param outfile:     the path to a file where the information in the celex dictionary is written in .json format.
    """

    with open(outfile, 'a+') as o_f:

        json.dump(celex_dict, o_f)


########################################################################################################################


def create_celex_dictionary(celex_dir, reduced=True, compounds=False):

    """
    :param celex_dir:   a string specifying the path to the folder where the information extracted from the Celex
                        database will be stored
    :param reduced:     a boolean specifying whether reduce phonological form should be always preferred when available
    :param compounds:   a boolean. If true, all entries in Celex are considered; if False, entries which contain spaces
                        are discarded
    :return celex_dict: a Python dictionary containing information about phonology and morphology about words from the
                        Celex database. For further details about which information is stored and how, check the
                        documentation of the functions in this module

    The function takes 6 seconds to run on a 2x Intel Xeon 6-Core E5-2603v3 with 2x6 cores and 2x128 Gb of RAM.
    """

    epw_path = os.path.join(celex_dir, 'epw.cd')
    epl_path = os.path.join(celex_dir, 'epl.cd')
    emw_path = os.path.join(celex_dir, 'emw.cd')
    eml_path = os.path.join(celex_dir, 'eml.cd')
    celex_vowels = vowels()
    celex_dict = initialize_celex_dict()

    celex_dict = read_epw(epw_path, celex_dict, celex_vowels, reduced=reduced, compounds=compounds)
    print(": ".join([strftime("%Y-%m-%d %H:%M:%S"), "I finished processing the epw.cd file."]))
    celex_dict = read_epl(epl_path, celex_dict)
    print(": ".join([strftime("%Y-%m-%d %H:%M:%S"), "I finished processing the epl.cd file."]))
    celex_dict = read_emw(emw_path, celex_dict, compounds=compounds)
    print(": ".join([strftime("%Y-%m-%d %H:%M:%S"), "I finished processing the emw.cd file."]))
    celex_dict = read_eml(eml_path, celex_dict)
    print(": ".join([strftime("%Y-%m-%d %H:%M:%S"), "I finished processing the eml.cd file."]))

    celex_dict = add_childes_words(celex_dict)

    dict_name = get_name(compounds=compounds, reduced=reduced)
    celex_dict_file = os.path.join(celex_dir, dict_name)

    print_celex_dict(celex_dict, celex_dict_file)

    return celex_dict


########################################################################################################################


def get_celex_dictionary(celex_dir, reduced=False, compounds=False):

    """
    :param celex_dir:       a string specifying the path to the Celex directory where the dictionary should be located.
                            If the dictionary is found, it is loaded and returned, otherwise it is created
    :param reduced:         a boolean specifying whether reduce phonological form should be always preferred when
                            available
    :param compounds:       a boolean. If true, all entries in Celex are considered; if False, entries which contain
                            spaces are discarded
    :return celex_dict:     a Python dictionary containing information about phonology and morphology of words extracted
                            from the CELEX database (for further details, see the documentation of the
                            celex_processing.py module
    """

    dict_name = get_name(compounds=compounds, reduced=reduced)

    try:
        celex_dict = json.load(open(os.path.join(celex_dir, dict_name), 'r'))
        print(": ".join([strftime("%Y-%m-%d %H:%M:%S"), "The desired Celex dictionary was loaded"]))
    except IOError:
        celex_dict = create_celex_dictionary(celex_dir, reduced=reduced, compounds=compounds)

    return celex_dict

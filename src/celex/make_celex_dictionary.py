__author__ = 'Gcassani'

"""This module creates a dictionary containing morphological and phonological information from the CELEX database,
   at both the token and type level, to be used to recode transcribed speech into phonologically and morphologically 
   richer representations (can be called from command line)"""

import os
import argparse
from src.celex.get import get_celex_dictionary


def main():

    parser = argparse.ArgumentParser(description='Process arguments to create Celex dictionary.')

    parser.add_argument('-C', '--Celex_dir', required=True, dest='celex_dir',
                        help='Specify the directory containing the CELEX files.')
    parser.add_argument('-r', '--reduced', action='store_true', dest='reduced',
                        help='Specify if the function considers reduced phonological forms.')

    args = parser.parse_args()

    if not os.path.isdir(args.celex_dir):
        raise ValueError("The folder you provided does not exist. Please, provide the path to an existing folder.")
    else:
        get_celex_dictionary(args.celex_dir, reduced=args.reduced)


########################################################################################################################


if __name__ == '__main__':

    main()

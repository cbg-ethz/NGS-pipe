# author = 'Simon Dirmeier'
# email  = 'simon.dirmeier@bsse.ethz.ch'
# date   = 07/01/16

from __future__ import print_function, absolute_import
import argparse
import subprocess


def parse_options():
    """
    Helper method to parse command line arguments. Returns the parsed command line arguments as string.

    :returns: returns the command line arguments as string
    :rtype: list(str)
    """
    parser = argparse.ArgumentParser(description='Parse a HMDB xml file.')
    parser.add_argument('-f', type=str, help='The data matrix you want to have analysed.', required=True,
                        metavar='input-file')
    parser.add_argument('-o', type=str, help='A path where all files are generated in.', required=True,
                        metavar='output-folder')
    parser.add_argument('-r', type=str, help='The path to the R-script', required=True,
                        metavar='output-folder')
    opts = parser.parse_args()
    return opts.r, opts.f, opts.o


def main():
    """
    Main method. Calls the Rscript! .... but maybe not, because calling R sucks!
    """

    # get parameters
    opts = parse_options()
    # initialize parser depending on chosen file format
    subprocess.call(['Rscript', opts[0], opts[1], opts[2]])


if __name__ == "__main__":
    main()

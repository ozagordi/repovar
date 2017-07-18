#!/usr/bin/env python3

def main():

    import sys
    import argparse
    from setuptools_scm import get_version

    from pkg_resources import get_distribution, DistributionNotFound
    try:
        __version__ = get_distribution('repovar').version
    except DistributionNotFound:
       # package is not installed
       pass


    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version',
                         version=__version__)
    # First options group
    group1 = parser.add_argument_group('Input files', 'Required input')
    group1.add_argument("-m", "--mutations", default="", type=str, dest="m",
                        help="csv file with annotated mutations")
    group1.add_argument("-s", "--subtype", default="", type=str, dest="s",
                        help="csv file with subtype evidence")

    args = parser.parse_args()

    # exit so that log file is not written
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    import logging
    import logging.handlers

    logging.basicConfig(filename='repovar.log', level=logging.DEBUG,
                        format='%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s', datefmt='%Y/%m/%d %H:%M:%S')
    logging.info(' '.join(sys.argv))

    from repovar import reportdrm
    reportdrm.main(mut_file=args.m, subtype_file=args.s)

if __name__ == "__main__": #  and __package__ is None:
    main()

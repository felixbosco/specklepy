import argparse


def parser(options=None):

    parser = argparse.ArgumentParser(description='This script searches for stars in the input file and creates a list file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--file', type=str, default=None, help='Fits file to search for stars.')
    parser.add_argument('-o', '--outfile', type=str, default=None, help='Name of the output file.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(options=None):

    args = parser(options=options)

    if args.file is None:
        logging.error("No file or file list was provided! Use --help for instructions.")
        raise RuntimeError("No file or file list was provided! Use --help for instructions.")

    


if __name__ == '__main__':
    main()

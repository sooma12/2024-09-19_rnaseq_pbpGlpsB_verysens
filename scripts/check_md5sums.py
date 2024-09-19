# check_md5sums.py
# MWS Jan 29th, 2024
"""
Take a SeqCenter xlsx file, e.g. "RNA Sequencing Stats.xlsx"
Also take a file containing file names and md5sum values for all files in the directory...
Print output that verifies md5sum values

*Run this on the cluster - Mac devices have md5 instead of md5sum, resulting in some wonkiness.
USAGE:
(first, run this to make the md5sum file: `md5sum * > gz_md5sums.txt`)
python check_md5sums.py <xlsx_file> <md5sum_file>

example:
python check_md5sums.py '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/input/Illumina RNA Reads/RNA Seq Statistics.xlsx' /work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/input/Illumina RNA Reads/gz_md5sums.txt

"""

import pandas as pd
import argparse

def get_args():
    """Return parsed CL arguments"""

    parser = argparse.ArgumentParser(
        description="Take SeqCenter xlsx file and md5sum text file and check values",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'xlsx_file',  # access with args.xlsx_file
        help="Provide filepath of xlsx input file",
        type=str)

    parser.add_argument(
        'md5sum_file',  #access with args.md5sum_file
        help="Provide filepath for output file; defaults to ./md5sum_file.txt",
        type=str)

    return(parser.parse_args())

def xlsx_md5sum_check(args):
    """ Read xlsx file and write files and their md5sums to new file"""
    df = pd.read_excel(args.xlsx_file)

    # Build list of md5 sum values
    md5sums_list = []
    for index, row in df.iterrows():
        md5sums_list.append(f"{row['R1 md5sum']}")
        md5sums_list.append(f"{row['R2 md5sum']}")

    # Read in file names and md5sum values from input file
    md5sums_pairs = []
    with open(args.md5sum_file, 'r') as md5fh:
        for line in md5fh:
            split_line = line.rstrip().split()
            split_line = [split_line[0], ' '.join(split_line[1:])]  # 2-member list = [md5sum, reassemble any files with spaces]
            pair = (split_line[0], split_line[1])
            md5sums_pairs.append(pair)

    # Check values
    for pair in md5sums_pairs:
        if pair[0] not in md5sums_list:
            print(f'No md5sum match for {pair[1]}.  Re-download this file.')
        else:
            print(f'{pair[1]}: OK')


if __name__ == "__main__":
    xlsx_md5sum_check(get_args())

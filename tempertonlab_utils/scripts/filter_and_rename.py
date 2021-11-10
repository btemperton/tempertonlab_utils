import argparse
import gzip

import numpy as np
from Bio import SeqIO


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--min_len", type=int, default=2000)
    parser.add_argument("--in_file", type=str, required=True)
    parser.add_argument("--out_file", type=str, required=True)
    parser.add_argument("--prefix", type=str)
    # Parse arguments
    args = parser.parse_args()
    return args

def get_N50(lengths):
    """
    Returns the N50 value for a list of lengths
    :param lengths: a list of lengths
    :return: The N50 value
    """
    sorted_list = np.sort(lengths)
    rev_sorted = sorted_list[::-1]
    mask = np.cumsum(rev_sorted) <= np.sum(lengths) / 2
    n50_id = np.sum(mask)
    return rev_sorted[n50_id]

def parse_file(handle):
    long = []
    old_lengths = []
    new_lengths = []
    counter = 1
    for record in SeqIO.parse(handle, 'fasta'):
        old_lengths.append(len(record.seq))
        if len(record.seq) >= args.min_len:
            record.description = f'old:{record.id}'
            record.name = ''
            new_lengths.append(len(record.seq))
            if args.prefix:
                record.id = f'{args.prefix}_{counter:06d}__len{len(record.seq)}bp'
                counter +=1
            long.append(record)
    return long, old_lengths, new_lengths

def main(args):

    try:
        with gzip.open(args.in_file, 'rt') as handle:
            long, old_lengths, new_lengths = parse_file(handle)
    except OSError:
        print(f'File {args.in_file} is not gzipped, processing as text')
        with open(args.in_file, 'r') as handle:
            long, old_lengths, new_lengths = parse_file(handle)

    if args.out_file.endswith('.gz'):
        with gzip.open(args.out_file, 'wt') as handle:
            SeqIO.write(long, handle, 'fasta')
    else:
        with open(args.out_file, 'w') as handle:
            SeqIO.write(long, handle, 'fasta')


    print(f'''{bcolors.OKBLUE}
Initial input file: {args.in_file}
Num contigs: {len(old_lengths)}
Min length: {np.min(old_lengths)}
Max length: {np.max(old_lengths)}
Median length: {np.median(old_lengths)}
N50: {get_N50(old_lengths)}
{bcolors.ENDC}
{bcolors.OKCYAN}
Filtered file: {args.out_file}
Num contigs: {len(new_lengths)}
Min length: {np.min(new_lengths)}
Max length: {np.max(new_lengths)}
Median length: {np.median(new_lengths)}
N50: {get_N50(new_lengths)}
{bcolors.ENDC}
''')


if __name__=="__main__":
    args = parse_args()
    main(args)
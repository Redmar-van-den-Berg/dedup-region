#!/usr/bin/env python3

import argparse
import logging


logging.basicConfig(
        format='%(message)s',
        level=logging.INFO,
)


def get_coverage(fname):
    """Return a list of coverage from fname"""
    cov = list()
    with open(fname) as fin:
        for line in fin:
            chrom, pos, depth = line.strip().split()
            cov.append(depth)
    return cov


def get_sample_name(fname):
    """Guess the sample name from the filename"""
    return fname.split('/')[0]


def main(args):
    # Read coverage files
    coverage = dict()
    for fname in args.input:
        name = get_sample_name(fname)
        coverage[name] = get_coverage(fname)

    # Print the coverage data
    samples = list(coverage.keys())

    # Determine the number of positions
    pos_count = len(coverage[samples[0]])

    # Print each sample
    print(*samples, sep='\t')
    for i in range(pos_count):
        print(*(coverage[s][i] for s in samples), sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('input', nargs='+')
    args = parser.parse_args()

    # Determine the output filenames
    main(args)

#!/usr/bin/env python3

import argparse
import logging
import os

import dnaio
import pysam
import xopen

logging.basicConfig(
        format='%(message)s',
        level=logging.INFO,
)


def get_readnames(bamfile):
    """ Get the readnames from the bamfile """
    readnames = set()

    with pysam.AlignmentFile(bamfile, 'rb') as samfile:
        for read in samfile.fetch():
            readnames.add(read.query_name)

    return readnames


def discard_from_fastq(read_names, fastqfile):
    """ Read the fastq file, discarding readnames as we go """
    fastq = dnaio.open(fastqfile, opener=xopen.xopen)
    for read in fastq:
        read_name = read.name.split(' ')[0]
        if read_name in read_names:
            read_names.remove(read.name)


def write_bam(infile, outfile, exclude):
    """ Write bame file, excluding the specified read names """
    with pysam.AlignmentFile(infile, 'rb') as bam_in:
        with pysam.AlignmentFile(outfile, 'wb', template=bam_in) as bam_out:
            for read in bam_in.fetch():
                if not read.query_name in exclude:
                    bam_out.write(read)


def main(args):
    # Read the bam file and save all read names
    read_names = get_readnames(args.input_bam)
    
    # Read the fastq file, discard all read_names that we encounter
    discard_from_fastq(read_names, args.input_fastq)

    # The reads that are left have been filtered out by umi-trie, so now we
    # output every read in bam_input, except those that are still present in
    # read_names
    write_bam(args.input_bam, args.output, read_names)


def output_file_name(infile, folder):
    """ Determine the output filename """
    fname = os.path.basename(infile)
    base, *rest = fname.split('.')
    base += '_dedup'
    fout = '.'.join([base] + rest)
    return f"{folder}/{fout}"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-fastq', required=True)
    parser.add_argument('--input-bam', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    # Determine the output filenames
    main(args)

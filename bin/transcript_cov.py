#!/usr/bin/env python3

import argparse
import dataclasses
import json

import gtf


@dataclasses.dataclass
class Depth:
    chrom: str
    pos: int
    depth: int


def get_transcript_positions(transcript):
    """Yield all positions that are part of the transcript"""
    positions = list()
    for exon in transcript:
        positions+=get_exon_positions(exon)
    return positions


def get_exon_positions(exon):
    """Yield all positions that are part of the exon"""
    # Get exon start and end as ints
    start = int(exon['start'])

    end = int(exon['end'])

    if exon['strand'] == '+':
        # We want to include the last bp of the exon
        for i in range(start, end + 1):
            yield i

    # Start position is before end on the chromosome, even for "-" strand
    # transcripts
    elif exon['strand'] == '-':
        # We want to include the 'last' bp of the exon, i.e. start-1
        for i in range(end, start - 1, -1):
            yield i


def get_transcript_depth(depths, transcript):
    """Yield the depth for every position within transcript"""
    # Make a dictionary of the depths so we can easily access the Depth objecs
    # by position
    d = {int(x.pos): x for x in depths}

    # Needed to create a dummy Depth object with 0 coverage
    chrom = depths[0].chrom

    for pos in get_transcript_positions(transcript):
        # Positions with 0 depth at the ends of the specified regions are not
        # included with samtools depth. Therefore, here we return 0 for missing
        # values
        yield d.get(pos, Depth(chrom, pos, 0))


def get_avg_exon_depth(depths, exon):
    """Calculate the average depth for each exon in transcript"""
    # The chromosome the exon is on
    chrom = exon['seqname']

    # The positions that are part of the exon
    positions = set(get_exon_positions(exon))

    # Get all depths that are part of the exon positions, on the correct
    # chromosome
    exon_depths = [x.depth for x in depths if x.chrom == chrom and x.pos in positions]

    return sum(exon_depths)/len(exon_depths)


def get_avg_transcript_depth(depths, transcript):
    return [get_avg_exon_depth(depths, exon) for exon in transcript]


def get_exons(gtf_fname, name, chrom_prefix):
    """ Extract the exons for transcript name from gtf """
    with open(gtf_fname) as fin:
        for record in gtf.gtf_to_json(fin, chrom_prefix, name):
            if record['feature'] == 'exon':
                if record['attribute']['transcript_id'] == name:
                    yield record


def get_coverage(cov_fname):
    """ Extract the coverage from the coverage file """
    coverage = list()
    with open(cov_fname) as fin:
        for line in fin:
            chrom, pos, depth = line.strip().split()
            coverage.append(Depth(chrom, int(pos), int(depth)))

    return coverage


def main(args):
    transcript = list(get_exons(args.gtf, args.transcript, 'chr'))
    coverage = get_coverage(args.coverage)

    # Print the transcript json
    if (out := args.transcript_out):
        with open(out, 'wt') as fout:
            print(json.dumps(transcript, indent=True), file=fout)

    # Print the coverage of the transcript
    if (out := args.coverage_out):
        with open(args.coverage_out, 'wt') as fout:
            for depth in get_transcript_depth(coverage, transcript):
                print(depth.chrom, depth.pos, depth.depth, sep='\t', file=fout)

    # Print the average exon coverage for each exon in transcript
    if (out := args.avg_exon):
        with open(out, 'wt') as fout:
            for depth in get_avg_transcript_depth(coverage, transcript):
                print(depth, file=fout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--gtf', required=True, help='GTF file')
    parser.add_argument('--transcript', required=True)
    parser.add_argument('--coverage', required=True)
    parser.add_argument('--coverage-out', required=False)
    parser.add_argument('--transcript-out', required=False)
    parser.add_argument('--avg-exon', required=False)

    arguments = parser.parse_args()
    main(arguments)

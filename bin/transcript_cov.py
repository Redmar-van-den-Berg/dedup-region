#!/usr/bin/env python3

import dataclasses

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
    d = {x.pos: x for x in depths}

    for pos in get_transcript_positions(transcript):
        yield d[pos]


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

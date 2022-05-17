#!/usr/bin/env python3

import dataclasses

@dataclasses.dataclass
class Depth:
    chrom: str
    pos: int
    depth: int


def get_transcript_positions(transcript):
    positions = list()
    for exon in transcript:
        positions+=get_exon_positions(exon)
    return positions


def get_exon_positions(exon):
    """ Yield all positions that are in the specified exon """
    # Get exon start and end as ints
    start = int(exon['start'])

    # We want to include the last bp of the exon
    end = int(exon['end'])

    if exon['strand'] == '+':
        for i in range(start, end + 1):
            yield i

    # Start position is before end on the chromosome, even for "-" strand
    # transcripts
    elif exon['strand'] == '-':
        for i in range(end, start - 1, -1):
            yield i

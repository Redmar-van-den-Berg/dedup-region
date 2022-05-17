#!/usr/bin/env python3

import dataclasses

@dataclasses.dataclass
class Depth:
    chrom: str
    pos: int
    depth: int


def get_transcript_positions(transcript):
    pass


def get_exon_positions(exon):
    """ Yield all positions that are in the specified exon """
    # Get exon start and end as ints
    start = int(exon['start'])

    # We want to include the last bp of the exon, so -1
    end = int(exon['end']) - 1

    if exon['strand'] == '+':
        for i in range(start, end + 2):
            yield i

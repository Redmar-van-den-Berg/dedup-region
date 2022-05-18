#!/usr/bin/env python3

import argparse
import gtf
import json

import dataclasses

@dataclasses.dataclass(order=True)
class Bed:
    chrom: str
    start: int
    end: int
    name: str

    def __repr__(self):
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}"


def transcript_to_bed(transcript):
    chrom = transcript['seqname']
    start = int(transcript['start']) - 1
    end = int(transcript['end']) - 1
    name = transcript['attribute']['transcript_id']
    return Bed(chrom, start, end, name)


def main(args):
    transcripts = set(args.transcripts)

    # Store the bed format, since we have no guarantee that the specified
    # transcripts are in order, but the bed output must be ordered
    bed = list()

    with open(args.gtf) as fin:
        for record in gtf.gtf_to_json(fin, chrom_prefix='chr'):
            # Stop if we found all transcripts
            if not transcripts:
                break

            if record['feature'] == 'transcript':
                if (transcript := record['attribute']['transcript_id']) in transcripts:
                    transcripts.remove(transcript)
                    bed.append(transcript_to_bed(record))
                
    for record in sorted(bed):
        print(record)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--transcripts', required=True, nargs='+')
    args = parser.parse_args()

    # Determine the output filenames
    main(args)

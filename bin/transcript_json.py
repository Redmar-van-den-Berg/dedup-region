#!/usr/bin/env python3

import argparse
import json

import gtf


def get_exons(gtf_fname, name, chrom_prefix):
    """ Extract the exons for transcript name from gtf """
    with open(gtf_fname) as fin:
        for record in gtf.gtf_to_json(fin, chrom_prefix, name):
            if record['feature'] == 'exon':
                if record['attribute']['transcript_id'] == name:
                    yield record


def main(args):
    transcript = list(get_exons(args.gtf, args.transcript, 'chr'))

    # Print the transcript json
    print(json.dumps(transcript, indent=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--gtf', required=True, help='GTF file')
    parser.add_argument('--transcript', required=True)

    arguments = parser.parse_args()
    main(arguments)

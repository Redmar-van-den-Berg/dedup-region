#!/usr/bin/env python3

import argparse
import plotly.graph_objects as go
import math

import gtf

def read_cov_depth(fname):
    depths = list()
    with open(fname) as fin:
        for line in fin:
            depths.append(float(line.strip()))
    return depths


def read_unique_reads(fname):
    """Read the number of unique reads from the umi-trie stats.dat file"""
    with open(fname) as fin:
        for line in fin:
            if line.startswith('unique:'):
                return int(line.split()[-1])


def guess_sample(fname):
    """ Guess the sample name from the file path """
    return fname.split('/')[0]


def guess_transcript(fname):
    """ Guess the transcript name from the file path """
    return fname.split('/')[-1].split('.')[0]


def guess_gene_name(fname, transcript):
    with open(fname) as fin:
        for record in gtf.gtf_to_json(fin, '', transcript):
            return record['attribute']['gene_name']


def main(args):
    # List of average exon depths
    data = list()

    # Corresponding list of sample names
    samples = list()

    # Corresponding list of unique reads
    unique_reads = list()

    for coverage, stats in zip(args.coverage, args.umi_stats):
        # Read the exon coverage from the filej
        cov = read_cov_depth(coverage)

        # Derive the sample name from the file name
        sample = guess_sample(stats)

        # Read the number of unique reads from the umi-tools stats file
        unique = read_unique_reads(stats)

        # Store to plot later
        data.append(cov)
        samples.append(sample)
        unique_reads.append(unique)


    # Guess the transcript name for the figure heading
    transcript_name = guess_transcript(coverage)

    # Get the gene name from the transcript name
    gene_name = guess_gene_name(args.gtf, transcript_name)

    display_name = f"{gene_name} ({transcript_name})"

    # Create the figure
    plot_coverage(data, samples, unique_reads, args.output, display_name,
            args.log_expression)


def get_max(data):
    """Get the max value from a nested list of integers"""
    m = 0
    for row in data:
        for value in row:
            m = max(m, value)
    return m


def plot_coverage(data , samples, unique_reads, fname, transcript_name, log_depth):
    # To use as x-axis
    exons = [f"Exon-{i+1}" for i in range(len(data[0]))]

    # Store the go.Scatter objects for every sample
    traces = list()

    for depths, sample, ureads in zip(data, samples, unique_reads):
        traces.append(go.Scatter3d(x=depths, y=exons, z=[ureads for _ in depths],
            name=sample, marker={'size': 1}))

    fig = go.Figure(
            data=traces
        )

    # Default orientation of the figure
    camera = dict(
            up=dict(x=1, y=0, z=0),
            eye=dict(x=0.5, y=0.5, z=-2)
    )

    # Determine the range for the x-axis plot
    max_value=get_max(data)
    if log_depth:
        # If all values are zero, we can't take log10
        if max_value == 0:
            x_range = [0, 1]
        # If not zero, we set that as the max value
        else:
            x_range = [0, math.log(max_value, 10)]
    else:
        x_range = [0, max_value]

    # Update the layout of the figure
    fig.update_layout(
            title=f"Exon coverage for {transcript_name}",
            scene_camera=camera,
            scene=dict(
                xaxis={
                    'title': 'Average exon coverage',
                    'type': 'log' if log_depth else 'linear',
                    'range': x_range
                },
                yaxis_title='',
                zaxis_title='Number of unique reads',
            )
    )
    fig.write_html(fname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--coverage', required=True, nargs='+', help='Coverage before deduplication')
    parser.add_argument('--umi-stats', required=True, nargs='+', help='umi-trie stats output file')
    parser.add_argument('--gtf', required=True, help='Used to determine the gene name')
    parser.add_argument('--output', required=True)
    parser.add_argument('--log-expression', default=False, action='store_true')

    arguments = parser.parse_args()
    main(arguments)

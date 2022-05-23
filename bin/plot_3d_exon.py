#!/usr/bin/env python3

import argparse
import plotly.graph_objects as go

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
    plot_coverage(data, samples, unique_reads, args.output, display_name)


def plot_coverage(data , samples, unique_reads, fname, transcript_name):
    print(data)
    print(samples)
    print(unique_reads)
    # To use as x-axis
    exons = [f"Exon-{i+1}" for i in range(len(data[0]))]

    # Store the go.Scatter objects for every sample
    traces = list()

    for data, sample, ureads in zip(data, samples, unique_reads):
        #traces.append(go.Scatter3d(y=data, x=exons, z=[ureads for _ in data], name=sample))
        traces.append(go.Scatter3d(x=data, y=exons, z=[ureads for _ in data], name=sample, marker={'size': 1}))
        #traces.append(go.Scatter3d(x=data, y=[ureads for _ in data], z=exons, name=sample, marker={'size': 1}))

    fig = go.Figure(
            data=traces
        )

    # Default orientation of the figure
    camera = dict(
            up=dict(x=1, y=0, z=0),
            eye=dict(x=1.5, y=1.5, z=-4)
    )
    fig.update_layout(
            title=f"Exon coverage for {transcript_name}",
            scene_camera=camera,
            scene=dict(
                xaxis_title='Average exon coverage',
                yaxis_title='',
                zaxis_title='Number of unique reads'
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

    arguments = parser.parse_args()
    main(arguments)

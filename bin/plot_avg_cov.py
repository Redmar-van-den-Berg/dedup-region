#!/usr/bin/env python3

import argparse
import plotly.graph_objects as go


def read_cov_depth(fname):
    depths = list()
    with open(fname) as fin:
        for line in fin:
            depths.append(float(line.strip()))
    return depths


def guess_sample(fname):
    """ Guess the sample name from the file path """
    return fname.split('/')[-1].split('.')[0]


def guess_transcript(fname):
    """ Guess the transcript name from the file path """
    return fname.split('/')[1]


def main(args):
    # List of percentage differences to plot
    diff = list()
    # Corresponding list of sample name
    samples = list() 
    for before, after in zip(args.before, args.after):
        # Some messing about to determine the difference in percentage between
        # before and after
        b = read_cov_depth(before)
        a = read_cov_depth(after)
        d = [ (af/be)* 100 for af, be in zip(a, b) ]

        # Derive the sample name from the file name
        sample = guess_sample(before)

        # Store to plot later
        diff.append(d)
        samples.append(sample)

    # Guess the transcript name for the figure heading
    transcript_name = guess_transcript(before)

    # Create the figure
    plot_coverage(diff, samples, args.output, transcript_name)


def plot_coverage(diff, samples, fname, transcript_name):
    # To use as x-axis
    exons = [f"Exon-{i+1}" for i in range(len(diff[0]))]

    # Store the go.Scatter objects for every sample
    traces = list()

    for sample, data in zip(samples, diff):
        traces.append(go.Scatter(y=data, x=exons, name=sample))

    fig = go.Figure(
            data=traces
        )

    fig.update_layout(
            title=f"Exon coverage for {transcript_name}",
            xaxis_title="Exon",
            yaxis_title="Coverage after deduplication(%)"
    )
    fig.write_html(fname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--before', required=True, nargs='+', help='Coverage before deduplication')
    parser.add_argument('--after', required=True, nargs='+', help='Coverage after deduplication')
    parser.add_argument('--output', required=True)

    arguments = parser.parse_args()
    main(arguments)

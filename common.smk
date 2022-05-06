containers = {
    "debian": "docker://debian:latest",
    # dnaio 0.8.1, pysam 0.19.0
    "dnaio": "docker://quay.io/biocontainers/mulled-v2-2996a7d035117c4238b2b801e740a69df21d91e1:6b3ae5f1a97f370227e8134ba3efc0e318b288c3-0",
    "samtools": "docker://quay.io/biocontainers/samtools:1.15.1--h1170115_0",
}

default = dict()


def get_fastq(wildcards, column):
    fastq = pep.sample_table.loc[wildcards.sample, column]

    # If a single fastq file is specified, forward will be a string
    if isinstance(fastq, str):
        return [fastq]
    # If multiple fastq files were specified, forward will be a list
    else:
        return fastq


def get_forward(wildcards):
    return get_fastq(wildcards, "forward")


def get_reverse(wildcards):
    return get_fastq(wildcards, "reverse")


def get_umi(wildcards):
    return get_fastq(wildcards, "umi")


def get_bamfile(wildcards):
    return pep.sample_table.loc[wildcards.sample, "bamfile"]

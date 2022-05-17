include: "common.smk"


pepfile: config["pepfile"]


# Apply the settings from the pepfile, overwriting the default ones
default.update(pep.config.get("dedup-region", dict()))

# Apply the options specified to snakemake, overwriting the default settings
# and the settings from the PEP file
default.update(config)

# Set the updated dict as the configuration for the pipeline
config = default


rule all:
    input:
        downsampled_bam=expand(
            "{sample}/regions/{sample}.dedup.bam",
            sample=pep.sample_table.sample_name,
        ),


rule concat:
    """Concatentate the input fastq files"""
    input:
        forw=get_forward,
        rev=get_reverse,
        umi=get_umi,
    output:
        forw="{sample}/concat/forward.fastq.gz",
        rev="{sample}/concat/reverse.fastq.gz",
        umi="{sample}/concat/umi.fastq.gz",
    log:
        "log/{sample}_concat.txt",
    container:
        containers["debian"]
    shell:
        """
        mkdir -p $(dirname {output.forw})

        cp {input.forw} {output.forw} || cat {input.forw} > {output.forw}
        cp {input.rev} {output.rev} || cat {input.rev} > {output.rev}
        cp {input.umi} {output.umi} || cat {input.umi} > {output.umi}
        """


rule umi_trie:
    """Run umi-trie on the fastq files"""
    input:
        forw=rules.concat.output.forw,
        rev=rules.concat.output.rev,
        umi=rules.concat.output.umi,
        umi_trie=config["umi_trie"],
    output:
        forw="{sample}/umi-trie/forward_dedup.fastq.gz",
        rev="{sample}/umi-trie/reverse_dedup.fastq.gz",
        umi="{sample}/umi-trie/umi_dedup.fastq.gz",
    log:
        "log/{sample}-umi-trie.txt",
    container:
        containers["dnaio"]
    shell:
        """
        folder=$(dirname {output.forw})
        mkdir -p $folder

        {input.umi_trie} \
            -d $folder \
            {input.forw} {input.rev} {input.umi} 2> {log}
        """


rule extract_regions:
    """Extract regions from the bam file"""
    input:
        bam=get_bamfile,
        bed=config["bedfile"],
    output:
        bam="{sample}/regions/{sample}.bam",
        bai="{sample}/regions/{sample}.bam.bai",
        cov="{sample}/regions/{sample}.cov",
    log:
        "log/{sample}/extract_regions.txt",
    container:
        containers["samtools"]
    shell:
        """
        samtools view \
            --bam \
            --with-header \
            --use-index \
            --regions-file {input.bed} \
            {input.bam} \
            > {output.bam} \
            2> {log}

            samtools index {output.bam} 2>> {log}
            samtools depth {output.bam} > {output.cov} 2>> {log}
        """


rule downsample_bam:
    """Downsample a bam file based on deduplicated fastq"""
    input:
        bam=rules.extract_regions.output.bam,
        fastq=rules.umi_trie.output.umi,
        downsample=srcdir("bin/downsample_bam.py"),
    output:
        bam="{sample}/regions/{sample}.dedup.bam",
        bai="{sample}/regions/{sample}.dedup.bam.bai",
        cov="{sample}/regions/{sample}.dedup.cov",
    log:
        "log/{sample}/downsample_bam.txt",
    container:
        containers["dnaio"]
    shell:
        """
        python3 {input.downsample} \
            --input-bam {input.bam} \
            --input-fastq {input.fastq} \
            --output {output.bam} 2> {log}

        # Index the output bam file
        python -c "import pysam; pysam.index('{output.bam}')" 2>> {log}
        python -c "import pysam; pysam.depth('-o', '{output.cov}', '{output.bam}')" 2>> {log}
        """

include: "common.smk"


pepfile: config["pepfile"]


# Apply the settings from the pepfile, overwriting the default ones
default.update(pep.config.get("dedup-region", dict()))

# Apply the options specified to snakemake, overwriting the default settings
# and the settings from the PEP file
default.update(config)

# Set the updated dict as the configuration for the pipeline
config = default

samples = pep.sample_table.sample_name


rule all:
    input:
        downsampled_bam=expand(
            "{sample}/regions/{sample}.dedup.bam",
            sample=samples,
        ),
        coverage=expand(
            "transcripts/{transcript}/{sample}.before.cov",
            sample=samples,
            transcript=config["transcripts"],
        ),


rule make_bedfile:
    """Make a bed file based on the specified transcripts"""
    input:
        gtf=config["gtf_file"],
        to_bed=srcdir("bin/transcript_to_bed.py"),
    output:
        bed="transcripts.bed",
    params:
        transcripts=config["transcripts"],
    log:
        "log/make_bedfile.txt",
    container:
        containers["dnaio"]
    shell:
        """
        python3 {input.to_bed} \
            --gtf {input.gtf} \
            --transcripts {params.transcripts} > {output.bed} 2> {log}
        """


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
        bed=rules.make_bedfile.output.bed,
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


rule extract_transcript_coverage:
    """Extract the coverage for each specified transcript"""
    input:
        gtf=config["gtf_file"],
        cov_before=rules.extract_regions.output.cov,
        cov_after=rules.downsample_bam.output.cov,
        transcript_cov=srcdir("bin/transcript_cov.py"),
    output:
        coverage_before="transcripts/{transcript}/{sample}.before.cov",
        exon_cov_before="transcripts/{transcript}/{sample}.before.avg.cov",
        coverage_after="transcripts/{transcript}/{sample}.after.cov",
        exon_cov_after="transcripts/{transcript}/{sample}.after.avg.cov",
    log:
        "log/{sample}/{transcript}/extract_transcript_coverage.txt",
    container:
        containers["dnaio"]
    shell:
        """
        # Create output folder
        mkdir -p transcripts/{wildcards.transcript}

        # Coverage before umi-tools
        python3 {input.transcript_cov} \
            --gtf {input.gtf} \
            --transcript {wildcards.transcript} \
            --coverage {input.cov_before} \
            --coverage-out {output.coverage_before} \
            --avg-exon {output.exon_cov_before} 2> {log}

        # Coverage after running umi-tools
        python3 {input.transcript_cov} \
            --gtf {input.gtf} \
            --transcript {wildcards.transcript} \
            --coverage {input.cov_after} \
            --coverage-out {output.coverage_after} \
            --avg-exon {output.exon_cov_after} 2>> {log}
        """

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


localrules:
    all,
    group_plots,


rule all:
    input:
        downsampled_bam=expand(
            "{sample}/regions/{sample}.dedup.bam",
            sample=samples,
        ),
        coverage=expand(
            "{sample}/{transcript}.before.cov",
            sample=samples,
            transcript=config["transcripts"],
        ),
        cov_plot=expand(
            "transcripts/{transcript}/{sample}_{transcript}.html",
            sample=samples,
            transcript=config["transcripts"],
        ),
        avg_plot=expand(
            "transcripts/{transcript}_exons.html",
            transcript=config["transcripts"],
        ),
        transcript_json=expand(
            "transcripts/{transcript}.json",
            transcript=config["transcripts"],
        ),
        transcript_after_merged=expand(
            "transcripts/{transcript}_after_merged.tsv",
            transcript=config["transcripts"],
        ),
        all_transcripts="transcripts/all_transcripts.html",
        all_3d_transcripts="transcripts/all_3d_transcripts.html",
        all_3d_transcripts_log="transcripts/all_3d_transcripts_log.html",


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
        containers["samtools"]
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
        stats="{sample}/umi-trie/stats.dat",
    log:
        "log/{sample}-umi-trie.txt",
    container:
        containers["dnaio"]
    shell:
        """
        folder=$(dirname {output.forw})
        mkdir -p $folder

        # umi-trie uses recursive functions
        ulimit -s 32768

        {input.umi_trie} \
            -d $folder \
            -s \
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


rule extract_transcript_json:
    """Extract the json for each specified transcript"""
    input:
        gtf=config["gtf_file"],
        transcript_json=srcdir("bin/transcript_json.py"),
    output:
        transcript_json="transcripts/{transcript}.json",
    log:
        "log/{transcript}/extract_transcript_coverage.txt",
    container:
        containers["dnaio"]
    shell:
        """
        # Create output folder
        mkdir -p transcripts/{wildcards.transcript}

        # Coverage before umi-tools
        python3 {input.transcript_json} \
            --gtf {input.gtf} \
            --transcript {wildcards.transcript} \
            > {output.transcript_json} 2> {log}
        """


rule extract_transcript_coverage:
    """Extract the coverage for each specified transcript"""
    input:
        gtf=config["gtf_file"],
        cov_before=rules.extract_regions.output.cov,
        cov_after=rules.downsample_bam.output.cov,
        transcript_cov=srcdir("bin/transcript_cov.py"),
    output:
        coverage_before="{sample}/{transcript}.before.cov",
        exon_cov_before="{sample}/{transcript}.before.avg.cov",
        coverage_after="{sample}/{transcript}.after.cov",
        exon_cov_after="{sample}/{transcript}.after.avg.cov",
    log:
        "log/{sample}/{transcript}/extract_transcript_coverage.txt",
    container:
        containers["dnaio"]
    shell:
        """
        # Create output folder
        mkdir -p {wildcards.sample}

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


rule plot_sample_coverage:
    """Plot the coverage, per exon, for each individual sample"""
    input:
        before=rules.extract_transcript_coverage.output.coverage_before,
        after=rules.extract_transcript_coverage.output.coverage_after,
        gtf=config["gtf_file"],
        plot_cov=srcdir("bin/plot_cov.py"),
    output:
        "transcripts/{transcript}/{sample}_{transcript}.html",
    log:
        "log/{sample}/plot_sample_coverage_{transcript}.txt",
    container:
        containers["plotly"]
    shell:
        """
        python3 {input.plot_cov} \
            --before {input.before} \
            --after {input.after} \
            --gtf {input.gtf} \
            --output {output} 2> {log}
        """


rule plot_average_transcript:
    """Plot the average difference between before and after for each sample"""
    input:
        before=get_average_coverage_before,
        after=get_average_coverage_after,
        gtf=config["gtf_file"],
        plot_avg_cov=srcdir("bin/plot_avg_cov.py"),
    output:
        "transcripts/{transcript}_exons.html",
    log:
        "log/plot_average_transcript_{transcript}.txt",
    container:
        containers["plotly"]
    shell:
        """
        python3 {input.plot_avg_cov} \
            --before {input.before} \
            --after {input.after} \
            --gtf {input.gtf} \
            --output {output} 2> {log}
        """


rule plot_3d_transcript:
    """Plot the 3d coverage per sample exon, against the number of unique reads"""
    input:
        coverage=get_average_coverage_after,
        stats=[f"{sample}/umi-trie/stats.dat" for sample in samples],
        gtf=config["gtf_file"],
        plot_3d_exon=srcdir("bin/plot_3d_exon.py"),
    output:
        linear="transcripts/{transcript}_3d.html",
        log="transcripts/{transcript}_3d_log.html",
    log:
        "log/plot_3d_transcript_{transcript}.txt",
    container:
        containers["plotly"]
    shell:
        """
        python3 {input.plot_3d_exon} \
            --coverage {input.coverage} \
            --umi-stats {input.stats} \
            --gtf {input.gtf} \
            --output {output.linear} 2> {log}

        python3 {input.plot_3d_exon} \
            --coverage {input.coverage} \
            --umi-stats {input.stats} \
            --gtf {input.gtf} \
            --log-expression \
            --output {output.log} 2>> {log}
        """


rule group_plots:
    """Join all average transcripts plots together"""
    input:
        exons=[f"transcripts/{t}_exons.html" for t in config["transcripts"]],
        exons_3d=[f"transcripts/{t}_3d.html" for t in config["transcripts"]],
        exons_3d_log=[f"transcripts/{t}_3d_log.html" for t in config["transcripts"]],
    output:
        exons="transcripts/all_transcripts.html",
        exons_3d="transcripts/all_3d_transcripts.html",
        exons_3d_log="transcripts/all_3d_transcripts_log.html",
    log:
        "log/group_plots.txt",
    container:
        containers["samtools"]
    shell:
        """
        cat {input.exons} > {output.exons} 2> {log}
        cat {input.exons_3d} > {output.exons_3d} 2>> {log}
        cat {input.exons_3d_log} > {output.exons_3d_log} 2>> {log}
        """


rule merge_to_tsv:
    """Merge the coverage for each exon per sample"""
    input:
        after=get_coverage_after,
        merge_cov=srcdir("bin/merge_transcript_cov.py"),
    output:
        "transcripts/{transcript}_after_merged.tsv",
    log:
        "log/merge_transcript_cov_{transcript}.txt",
    container:
        containers["plotly"]
    shell:
        """
        python3 {input.merge_cov} \
            {input.after} > {output} 2> {log}
        """

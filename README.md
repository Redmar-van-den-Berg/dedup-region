[![Continuous Integration](https://github.com/Redmar-van-den-Berg/dedup-region/actions/workflows/ci.yml/badge.svg)](https://github.com/Redmar-van-den-Berg/dedup-region/actions/workflows/ci.yml)
[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)
![GitHub release](https://img.shields.io/github/v/release/redmar-van-den-berg/dedup-region)
![Commits since latest release](https://img.shields.io/github/commits-since/redmar-van-den-berg/dedup-region/latest)

# dedup-region
Example of a snakemake project

## Installation
Download the repository from github
```bash
git clone https://github.com/Redmar-van-den-Berg/dedup-region.git
```

Install and activate the
[conda](https://docs.conda.io/en/latest/miniconda.html)
environment.
```bash
conda env create --file environment.yml
conda activate dedup-region
```

## Settings
There are three levels where configuration options are set, in decreasing order
of priority.
1. Flags passed to snakemake using `--config`, or in the specified
   `--configfile`.
2. Setting specified in the PEP project configuration, under the key
   `dedup-region`
3. The default settings for the pipeline, as specified in the `common.smk` file

### Supported settings
The following settings are available for the pipeline.
| Option               | Type              | Explanation                             |
| ---------------------| ----------------- | --------------------------------------- |
| umi_trie             | Required binary   | Compiled version of umi-trie (the tests use a dummy version that is not suitable for real data |
| gtf file             | Required gtf file | http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz |
| transcripts          | List of transcripts | Must match the `transcript_id` from the gtf file |

## Tests
You can run the tests that accompany this pipeline with the following commands

```bash
# Check if requirements are installed, and run linting on the Snakefile
pytest --kwd --tag sanity

# Test the pipeline settings in dry-run mode
pytest --kwd --tag dry-run

# Test the performance of the pipeline by running on the test data
pytest --kwd --tag integration
```

## Limitations
* For transcripts located on the mitochondrion, make sure that you rename the
entries in the GTF file to use `M` instead of `MT` as chromosome name.
* If a transcript has no coverage **before** deduplication, the remaining coverage after deduplication will be shown as 0%.
0% in the `trasncripts/{transcript}_exons.html graph.

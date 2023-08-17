
from math import factorial

#snakemake --cores 28  kmer_run --latency-wait 60 --rerun-incomplete -p --verbose --configfile=config/config_paper.yaml
#snakemake --cores 28  pantools_run --latency-wait 60 --rerun-incomplete -p --verbose --configfile=config/config_paper.yaml
#snakemake --cores 28  pantools_panmatrix_run --latency-wait 60 --rerun-incomplete -p --verbose --configfile=config/config_paper.yaml
#snakemake --cores 28  roary_run --latency-wait 60 --rerun-incomplete -p --verbose --configfile=config/config_paper.yaml

#--------------------------------------------------------------------------------
# INPUT
#--------------------------------------------------------------------------------
OUTDIR=config["OUTDIR"]
DATA=config["DATA"]

#--------------------------------------------------------------------------------
# RULES
#--------------------------------------------------------------------------------
include: "rules/util.smk"
include: "rules/fasta_stats.smk"
#--------------------------------------------------------------------------------
include: "rules/kmer.smk"
include: "rules/prokka.smk"
include: "rules/roary.smk"
include: "rules/pantools.smk"
#--------------------------------------------------------------------------------

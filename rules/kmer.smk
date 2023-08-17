rule kmer_run:
    input:
        expand(OUTDIR+"data/kmer_{k}/{species}/pangenome_K:{k}_can.txt", species=config["pangenome"],k=config["K"]),

rule kmer_plot:
    input:
        expand(OUTDIR+"plots/{species}/new_kmer_opt_k:{k}.png", species=config["pangenome"], k=config["K"]),
        expand(OUTDIR+"plots/{species}/new_kmer_opt_log_k:{k}.png", species=config["pangenome"], k=config["K"]),
        expand(OUTDIR+"plots/{species}/total_kmer_opt_k:{k}.png", species=config["pangenome"], k=config["K"]),
        expand(OUTDIR+"plots/{species}/total_kmer_opt_log_k:{k}.png",species=config["pangenome"], k=config["K"])
        
#--------------------------------------------------------------------------------

## k-mer count
# Automatic behavior: count kmers in canonincal form. 
# Use -b to deactivate
rule kmer_count_hist:
    input:
        fasta_path
    output:
        OUTDIR+"data/kmer_{k}/{species}/hist_{k}_can.txt"
    threads:
        10
    shell:
        "./scripts/pangrowth/pangrowth hist -k{wildcards.k} -t{threads} {input} > {output}"

rule kmer_plot_out:
    input:
        OUTDIR+"data/kmer_{k}/{species}/pangenome_K:{k}_can.txt"
    output:
        OUTDIR+"plots/{species}/new_kmer_opt_k:{k}.png", 
        OUTDIR+"plots/{species}/new_kmer_opt_log_k:{k}.png", 
        OUTDIR+"plots/{species}/total_kmer_opt_k:{k}.png", 
        OUTDIR+"plots/{species}/total_kmer_opt_log_k:{k}.png"
    shell:
        "Rscript scripts/create_kmer_plots.R {wildcards.species} {wildcards.k}"

## PANGENOME SIZE
rule pangenome_size_kmer:
    input:
        OUTDIR+"data/kmer_{k}/{species}/hist_{k}_can.txt"
    output:
        OUTDIR+"data/kmer_{k}/{species}/pangenome_K:{k}_can.txt"
    shell:
        "./scripts/pangrowth/pangrowth growth -h {input} > {output}"

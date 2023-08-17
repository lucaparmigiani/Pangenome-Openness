GENE_ROARY_DATA=OUTDIR+"data/gene_roary"
GENE_ROARY_PLOT=OUTDIR+"plot/gene_roary"

rule roary_run:
    input:
        expand(GENE_ROARY_DATA+"/{species}/pangenome_gene.txt", species=config["pangenome"]),

rule roary_plot:
    input:
        expand(GENE_ROARY_PLOT+"/{species}/new_gene_opt.png", species=config["pangenome"]),
        expand(GENE_ROARY_PLOT+"/{species}/new_gene_opt_log.png", species=config["pangenome"]),
        expand(GENE_ROARY_PLOT+"/{species}/total_gene_opt.png", species=config["pangenome"]),
        expand(GENE_ROARY_PLOT+"/{species}/total_gene_opt_log.png",species=config["pangenome"]),
#--------------------------------------------------------------------------------

rule roary_start:
    input:
        gff_path
        #expand("prokka/{species}/{sample}", sample=[1,23])
    output:
        OUTDIR+"roary/{species}/number_of_genes_in_pan_genome.Rtab",
        OUTDIR+"roary/{species}/number_of_new_genes.Rtab",
        OUTDIR+"roary/{species}/gene_presence_absence.Rtab"
    benchmark:
        "benchmarks/{species}/roary/roary_run.benchmark.txt"
    threads:
        28
    shell:
        #Split paralogous (default)
        "rm -r "+OUTDIR+"roary/{wildcards.species}/; roary -p {threads} -f "+OUTDIR+"roary/{wildcards.species} -n -v {input}"
        #Do not split paralogous (-s)
        #"rm -r "+OUTDIR+"roary/{wildcards.species}/; roary -s -p {threads} -f "+OUTDIR+"roary/{wildcards.species} -n -v {input}"

rule roary_move:
    input:
        x=OUTDIR+"roary/{species}/number_of_genes_in_pan_genome.Rtab",
        y=OUTDIR+"roary/{species}/number_of_new_genes.Rtab",
        z=OUTDIR+"roary/{species}/gene_presence_absence.Rtab"
    output:
        x=GENE_ROARY_DATA+"/{species}/number_of_genes_in_pan_genome.Rtab",
        y=GENE_ROARY_DATA+"/{species}/number_of_new_genes.Rtab",
        z=GENE_ROARY_DATA+"/{species}/gene_presence_absence.Rtab"
    shell:
        "cp {input.x} {output.x}; cp {input.y} {output.y}; cp {input.z} {output.z}"

rule clean_panmatrix:
    input:
       GENE_ROARY_DATA+"/{species}/gene_presence_absence.Rtab"
    output:
       GENE_ROARY_DATA+"/{species}/gene_presence_absence_clean.Rtab"
    shell:
        """
        tail +2 {input} | cut -f2- > {output}
        sed -i 's/\t/ /g' {output} 
        """

rule roary_plots:
    input:
        GENE_ROARY_DATA+"/{species}/number_of_genes_in_pan_genome.Rtab",
        GENE_ROARY_DATA+"/{species}/number_of_new_genes.Rtab",
        GENE_ROARY_DATA+"/{species}/pangenome_gene.txt"
    output:
        GENE_ROARY_PLOT+"/{species}/total_gene.png",
        GENE_ROARY_PLOT+"/{species}/new_gene.png",
        GENE_ROARY_PLOT+"/{species}/new_gene_opt.png",
        GENE_ROARY_PLOT+"/{species}/new_gene_opt_log.png",
        GENE_ROARY_PLOT+"/{species}/total_gene_opt.png", 
        GENE_ROARY_PLOT+"/{species}/total_gene_opt_log.png"
    shell:
        "Rscript scripts/create_gene_plots.R {wildcards.species}"

rule pangenome_size_gene:
    input:
        GENE_ROARY_DATA+"/{species}/gene_presence_absence_clean.Rtab"
    output:
        GENE_ROARY_DATA+"/{species}/pangenome_gene.txt"
    shell:
        "./scripts/pangrowth/pangrowth growth -p {input} > {output}"

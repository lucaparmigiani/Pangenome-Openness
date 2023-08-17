GENE_PANTOOLS_DATA=OUTDIR+"data/gene_pantools"
GENE_PANTOOLS_PLOT=OUTDIR+"plot/gene_pantools"

rule pantools_run:
    input:
        expand(OUTDIR+"pantools/db/{species}/done.txt", species=config["pangenome"]),

rule pantools_panmatrix_run:
    input:
        expand(OUTDIR+"pantools/db/{species}/gene_presence_absence.txt", species=config["pangenome"]),
#--------------------------------------------------------------------------------

rule pantools_GFF:
    input:
        gff_path
    output:
        OUTDIR+"pantools/gff/{species}_GFF.txt"
    run:
        with open(str(output), "w") as out:
            i = 1
            for gff in input:
                out.write(str(i)+" "+gff+"\n")
                i+=1

rule pantools_FA:
    input:
        fasta_path
    output:
        OUTDIR+"pantools/fa/{species}_fa.txt"
    run:
        with open(str(output), "w") as out:
            for fa in input:
                out.write(fa+"\n")

rule pantools_start:
    input: 
        gff=OUTDIR+"pantools/gff/{species}_GFF.txt",
        fa=OUTDIR+"pantools/fa/{species}_fa.txt", 
    output:
        OUTDIR+"pantools/db/{species}/done.txt"
    params:
        #pantools="java -Xms50g -Xmx50g -jar scripts/pantools/target/pantools-3.4.jar",
        pantools="java -Xms50g -Xmx50g -jar /home/ubuntu/Tools/pantools/pantools-4.2.2.jar",
        dir_species=OUTDIR+"pantools/db/{species}",
        log=OUTDIR+"pantools/log/{species}.txt",
        logdir=OUTDIR+"pantools/log",
        relaxation=1,
    benchmark:
        "benchmarks/{species}/pantools/pantools_run.benchmark.txt"
    threads:
        28
    shell:
        """
        mkdir -p {params.logdir}
        echo 'started' > {params.log}
        echo '{wildcards.species}' | tee -a {params.log}
        echo 'build_pangenome' | tee -a {params.log}
        {params.pantools} build_pangenome {params.dir_species} {input.fa} -t {threads}
        echo 'add_annotations' | tee -a {params.log}
        {params.pantools} add_annotations {params.dir_species} {input.gff}
        echo 'group' | tee -a {params.log}
        {params.pantools} group {params.dir_species} -t {threads} --relaxation {params.relaxation}
        echo 'pangenome_structure (genes)' | tee -a log.txt
        {params.pantools} pangenome_structure --loops 500 -t {threads} {params.dir_species}
        cd {params.dir_species}/pangenome_size/gene
        Rscript heaps_law.R
        cd ../../../../../..
        echo '1' > {output}
        """

rule pantools_panmatrix:
    input: 
        OUTDIR+"pantools/db/{species}"
    output:
        OUTDIR+"pantools/db/{species}/gene_presence_absence.txt"
    params:
        #pantools="java -Xms50g -Xmx50g -jar scripts/pantools/target/pantools-3.4.jar",
        pantools="java -Xms50g -Xmx50g -jar /home/ubuntu/Tools/pantools/pantools-4.2.2.jar",
        dir_species=OUTDIR+"pantools/db/{species}",
        raw_panmatrix=OUTDIR+"pantools/db/{species}/gene_classification/classified_groups.csv"
    shell:
        """
        {params.pantools} gene_classification {params.dir_species}
        bash scripts/pantools_panmatrix.sh {params.raw_panmatrix} > {output}
        """
        
#rule pantools_move:
#    input:
#        x=OUTDIR+"pantools/{species}/number_of_genes_in_pan_genome.Rtab",
#        y=OUTDIR+"pantools/{species}/number_of_new_genes.Rtab",
#        z=OUTDIR+"pantools/{species}/gene_presence_absence.Rtab"
#    output:
#        x=GENE_PANTOOLS_DATA+"/{species}/number_of_genes_in_pan_genome.Rtab",
#        y=GENE_PANTOOLS_DATA+"/{species}/number_of_new_genes.Rtab",
#        z=GENE_PANTOOLS_DATA+"/{species}/gene_presence_absence.Rtab"
#    shell:
#        "cp {input.x} {output.x}; cp {input.y} {output.y}; cp {input.z} {output.z}"
#
#rule clean_panmatrix:
#    input:
#       GENE_PANTOOLS_DATA+"/{species}/gene_presence_absence.Rtab"
#    output:
#       GENE_PANTOOLS_DATA+"/{species}/gene_presence_absence_clean.Rtab"
#    shell:
#        """
#        tail +2 {input} | cut -f2- > {output}
#        sed -i 's/\t/ /g' {output} 
#        """
#
##rule pantools_plots:
##    input:
##        GENE_PANTOOLS_DATA+"/{species}/number_of_genes_in_pan_genome.Rtab",
##        GENE_PANTOOLS_DATA+"/{species}/number_of_new_genes.Rtab",
##        GENE_PANTOOLS_DATA+"/{species}/pangenome_gene.txt"
##    output:
##        GENE_PANTOOLS_PLOT+"/{species}/total_gene.png",
##        GENE_PANTOOLS_PLOT+"/{species}/new_gene.png",
##        GENE_PANTOOLS_PLOT+"/{species}/new_gene_opt.png",
##        GENE_PANTOOLS_PLOT+"/{species}/new_gene_opt_log.png",
##        GENE_PANTOOLS_PLOT+"/{species}/total_gene_opt.png", 
##        GENE_PANTOOLS_PLOT+"/{species}/total_gene_opt_log.png"
##    shell:
##        "Rscript scripts/create_gene_plots.R {wildcards.species}"
#
#rule pangenome_size_gene:
#    input:
#        GENE_PANTOOLS_DATA+"/{species}/gene_presence_absence_clean.Rtab"
#    output:
#        GENE_PANTOOLS_DATA+"/{species}/pangenome_gene.txt"
#    shell:
#        "./scripts/pangrowth/pangrowth growth -p -i {input} > {output}"


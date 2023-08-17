rule prokka_annotation:
    input:
        fna=DATA+"{species}/{sample}.fna"
    output:
        directory(OUTDIR+"prokka/{species}/{sample}")
    params:
        genus=get_genus,
        kingdom=get_kingdom
    threads:
        14 
    benchmark:
        "benchmarks/{species}/prokka/prokka_annotation_{sample}.benchmark.txt"
    #conda: "envs/prokka.yaml"
    #container: "docker://staphb/prokka:latest"
    # --singularity-args "-B $PWD:/data"
    shell:
        "prokka --kingdom {params.kingdom} "
        "--outdir "+OUTDIR+"/prokka/{wildcards.species}/{wildcards.sample} "
        "--genus {params.genus}  "
        "{input.fna} --cpus {threads} --force"

rule prokka_run_and_clean:
    input:
        OUTDIR+"prokka/{species}/{sample}"
    output:
        OUTDIR+"prokka/{species}/{sample}.gff"
    shell:
        "cp {input}/PROKKA_[0-9]*.gff {output}"

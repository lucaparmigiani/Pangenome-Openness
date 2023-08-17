suppressWarnings(library(dplyr))
suppressWarnings(library(tidyverse))
source("soft_processing.R")

#--------------------------------------------------------------------------------
#DATA FOLDERS
folder = "data"

folder_seqkit = file.path(folder,"fasta_stats")

kmer_path = file.path(folder, "kmer")

pantools_path = file.path(folder, "pantools", "relaxation")

roary_path = file.path(folder, "roary")
roary_no_split_path = file.path(folder, "roary_no_split")

prokka_path = file.path(folder, "panmatrix/prokka/growth")

bpga_path = file.path(folder, "bpga")

folder_benchmarks = file.path(folder, "benchmarks")

# Other folders
mat_pantools_path = file.path(folder, "panmatrix/pantools/growth_relaxation")
mat_pantools_hist_path = file.path(folder, "panmatrix/pantools/panmatrix_relaxation")
mat_roary_path = file.path(folder, "panmatrix/roary/growth")
mat_roary_hist_path = file.path(folder, "panmatrix/roary/panmatrix")

#--------------------------------------------------------------------------------
# GENERAL
#--------------------------------------------------------------------------------
#pantools since it has the least amount of species we are sure we have data for
#every other type
species <- function(type="pantools") {
    species = c()
    if (type == "roary") {
        species = list.files(roary_path)
    } 
    else if (type == "pantools") {
        species = list.files(paste0(pantools_path,1)) #1 for relaxation 1
    } 
    else if (type == "kmer"){
        species = list.files(kmer_path)
    } 
    return(species)
}

species_df <- function(type="pantools") {
    species = c()
    if (type == "roary") {
        species = list.files(roary_path)
    } 
    else if (type == "pantools") {
        species = list.files(paste0(pantools_path,1)) #1 for relaxation 1
    } 
    else if (type == "kmer") {
        species = list.files(kmer_path)
    }

    df = as.data.frame(matrix(nrow=0,ncol=0))
    for (sp in species) {
        df = rbind(df, c(sp, unlist(summarise_seqkit(file.path(folder_seqkit, sp, "fasta_stats.csv")))))
    }
    colnames(df) = c("species",colnames(summarise_seqkit(file.path(folder_seqkit, species[1], "fasta_stats.csv"))))
    df[,2:ncol(df)] = df[,2:ncol(df)] %>% mutate_all(as.numeric)

    return(df)
}

species_all_items = function() {
    kmer = get_num_items("kmer")
    roary = get_num_items("roary")
    pantools = get_num_items("pantools")
    bpga = get_num_items("bpga")
    df = list(kmer, roary, pantools, bpga) %>% reduce(full_join, by="species")
    colnames(df)=c( "Tot.Pangrowth","species", "Single.Pangrowth",
                    "Tot.Roary","Single.Roary",
                    "Tot.Pantools","Single.Pantools",
                    "Tot.BPGA", "Single.BPGA")

    return(species_df() %>% full_join(df, by="species"))
}

read_avg_new = function (sp, program, ...) {
    ret = NA
    if (program=="kmer") {
        ret = read_kmer(sp, ...)
    }
    else if (program=="roary") {
        ret = read_roary(sp, ...)
    }
    else if (program=="roary_mat") {
        ret = read_mat_roary(sp, ...)
    }
    else if (program=="pantools") {
        ret = read_pantools(sp, ...)
    }
    else if (program=="pantools_mat") {
        ret = read_mat_pantools(sp, ...)
    }
    else if (program=="bpga") {
        ret = read_bpga_new(sp, ...)
    }
    else if (program=="bpga_mat") {
        ret = read_mat_bpga_new(sp, ...)
    } else{
        cat("Unknown program\n")
    }
    return(ret)
}

read_avg_tot = function (sp, program, ...) {
    ret = NA
    if (program=="kmer") {
        ret = read_kmer_tot(sp, ...)
    }
    #else if (program=="roary") {
    #    ret = read_roary_tot(sp, ...)
    #}
    #else if (program=="pantools") {
    #    ret = read_pantools(sp, ...)
    #}
    #else if (program=="bpga") {
    #    ret = read_bpga_new(sp, ...)
    #}
    else if (program=="roary_mat") {
        ret = read_mat_roary_tot(sp, ...)
    }
    else if (program=="pantools_mat") {
        ret = read_mat_pantools_tot(sp, ...)
    }
    else if (program=="bpga_mat") {
        ret = read_mat_bpga_tot(sp, ...)
    } else{
        cat("Unknown program\n")
    }
    return(ret)
}

#--------------------------------------------------------------------------------
# INFOS
#--------------------------------------------------------------------------------
get_num_items <- function (program, ...) {  
    items = c()
    bacteria = c()
    single_genome = c()
    for (sp in species()) {
        if (program == "roary_mat" | program == "roary") {
            res = read_roary_hist(sp, ...)
            sng = read_mat_roary_tot(sp, ...)$y[1]
        } 
        else if (program == "pantools_mat" | program == "pantools") {
            res = read_pantools_hist(sp,...)
            sng = read_mat_pantools_tot(sp, ...)$y[1]
        } 
        else if (program == "kmer"){
            res = read_kmer_hist(sp,...)
            sng = read_kmer_tot(sp, ...)$y[1]
        } 
        else if (program == "bpga_mat" | program == "bpga"){
            res = read_bpga_hist(sp, ...)
            sng = read_mat_bpga_tot(sp, ...)$y[1]
        }
        num = sum(res)
        items = c(items, num)
        bacteria = c(bacteria, sp)
        single_genome = c(single_genome, sng)
    }
    return(data.frame(items=items, species=bacteria, single_genome=single_genome))
}


#--------------------------------------------------------------------------------
# Seqkit
#--------------------------------------------------------------------------------
read_seqkit = function (filepath) {
    read.csv(filepath)
}

summarise_seqkit = function(filepath) {
    df = read_seqkit(filepath)
    #head(df)
    return(data.frame(
        N = nrow(df),
        avg_sum_len = mean(df$sum_len),
        avg_N50 = mean(df$N50),
        tot_bp = sum(df$sum_len),
        min_len = min(df$min_len),
        max_len = min(df$max_len)
    ))
}


#--------------------------------------------------------------------------------
# PANTOOLS
#--------------------------------------------------------------------------------
read_pantools <- function (sp, relaxation=1, avg=TRUE) {
    df1 <- read.csv(file.path(paste0(pantools_path,relaxation),sp,
                              "pangenome_size/gene/gene_nrs_for_heaps_law.csv"), header=FALSE)
    df2 <- read.csv(file.path(paste0(pantools_path,relaxation), sp, 
                              "pangenome_size/gene/genes_for_heaps_law.csv"), header = FALSE)
    x <- as.integer(df1)
    y <- as.numeric(df2)

    # Remove x=1
    y <- y[x!=1]
    x <- x[x!=1]

    if (avg) {
        Y = matrix(y, ncol=(max(x)-min(x)+1), byrow=T)
        Y = colMeans(Y)
        X = min(x):max(x)
    } else {
        X = x
        Y = y
    }

    p0 <- c(mean(Y[which(X == 2)] ), 1)
    return(list(x=X,y=Y,p0=p0))
}

read_pantools_tot <- function (sp, relaxation=1) {
    df1 <- read.csv(file.path(paste0(pantools_path, relaxation),sp,
                              "pangenome_size/gene/gene_nrs_for_heaps_law.csv"), header=FALSE)
    df2 <- read.csv(file.path(paste0(pantools_path, relaxation), sp, 
                              "pangenome_size/gene/genes_for_heaps_law.csv"), header = FALSE)
    x <- as.integer(df1)
    y <- as.numeric(df2)

    pan = read.csv(file.path(paste0(pantools_path, relaxation),sp,"gene_classification/classified_groups.csv"))
    pan[ncol(pan)] = NULL
    panmatrix = as.matrix(pan[,3:ncol(pan)])

    y1 = as.vector(colSums(panmatrix))
    x1 = rep(1, length(y1))

    x = c(x1,x)
    y = c(y1,y)

    p0 <- c(mean(y[which(x == 1)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_mat_pantools <- function (sp, relaxation=1) {
    input_path = file.path(paste0(mat_pantools_path,relaxation), sp)
    tot_gene = as.vector(unlist(read.csv(input_path, header=F, sep=" ")))
    tot_gene = tot_gene[!is.na(tot_gene)]
    new_gene = vector("numeric", length(tot_gene))
    new_gene[1] = tot_gene[1]
    for (i in 2:(length(tot_gene))) {
        new_gene[i] = tot_gene[i] - tot_gene[i-1]
    }

    N = length(new_gene)
    df2 = new_gene
    df1 = 1:N

    x <- as.integer(df1)
    y <- as.numeric(df2)
    # Remove x=1
    y <- y[x!=1]
    x <- x[x!=1]

    p0 <- c(mean(y[which(x == 2)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_mat_pantools_tot <- function (sp, relaxation=1) {
    input_path = file.path(paste0(mat_pantools_path,relaxation), sp)
    tot_gene = as.vector(unlist(read.csv(input_path, header=F, sep=" ")))
    tot_gene = tot_gene[!is.na(tot_gene)]

    N = length(tot_gene)
    df2 = tot_gene
    df1 = 1:N

    x <- as.integer(df1)
    y <- as.numeric(df2)

    p0 <- c(mean(y[which(x == 1)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_pantools_hist <- function (sp, relaxation=1) {
    input_path =file.path(paste0(mat_pantools_hist_path,relaxation), sp, "panmatrix.txt")
    panmatrix = read.csv(input_path, header=F, sep=" ")
    counts = table(rowSums(panmatrix))
    N = ncol(panmatrix)
    tab <- rep(0, N)
    names(tab) <- 1:N
    tab[names(counts)] <- counts
    tab
}
#--------------------------------------------------------------------------------
# ROARY
#--------------------------------------------------------------------------------
read_roary <- function (sp, avg=TRUE) {
    df2 <- read.table(file.path(roary_path, sp, "number_of_new_genes.Rtab"), header = FALSE, sep="\t")
    N = ncol(df2)
    perm = nrow(df2)


    df2 = c(t(df2))
    df1 <- rep(c(1:N), perm)
    x <- as.integer(df1)
    y <- as.numeric(df2)
    # Remove x=1
    y <- y[x!=1]
    x <- x[x!=1]

    if (avg) {
        Y = matrix(y, ncol=(max(x)-min(x)+1), byrow=T)
        Y = colMeans(Y)
        #Y = matrixStats::colMedians(Y)
        X = min(x):max(x)
    } else {
        X = x
        Y = y
    }

    p0 <- c(mean(Y[which(X == 2)] ), 1)
    return(list(x=X,y=Y,p0=p0))
}
read_roary_no_split <- function (sp) {
    df2 <- read.table(file.path(roary_no_split_path, sp, "number_of_new_genes.Rtab"), header = FALSE, sep="\t")
    N = ncol(df2)
    perm = nrow(df2)
    df2 = c(t(df2))
    df1 <- rep(c(1:N), perm)
    x <- as.integer(df1)
    y <- as.numeric(df2)
    # Remove x=1
    y <- y[x!=1]
    x <- x[x!=1]

    p0 <- c(mean(y[which(x == 2)] ), 1)
    return(list(x=x,y=y,p0=p0))
}
read_roary_tot <- function (sp) {
    df2 <- read.table(file.path(roary_path, sp, "number_of_genes_in_pan_genome.Rtab"), header = FALSE, sep="\t")
    N = ncol(df2)
    perm = nrow(df2)
    df2 = c(t(df2))
    df1 <- rep(c(1:N), perm)
    x <- as.integer(df1)
    y <- as.numeric(df2)
    head(y)

    p0 <- c(mean(y[which(x == 1)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_mat_roary <- function (sp) {
    input_path =file.path(mat_roary_path, sp)
    tot_gene = as.vector(unlist(read.csv(input_path, header=F, sep=" ")))
    tot_gene = tot_gene[!is.na(tot_gene)]
    new_gene = vector("numeric", length(tot_gene))
    new_gene[1] = tot_gene[1]
    for (i in 2:(length(tot_gene))) {
        new_gene[i] = tot_gene[i] - tot_gene[i-1]
    }

    N = length(new_gene)
    df2 = new_gene
    df1 = 1:N

    x <- as.integer(df1)
    y <- as.numeric(df2)
    # Remove x=1
    y <- y[x!=1]
    x <- x[x!=1]

    p0 <- c(mean(y[which(x == 2)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_mat_roary_tot <- function (sp) {
    input_path =file.path(mat_roary_path, sp)
    tot_gene = as.vector(unlist(read.csv(input_path, header=F, sep=" ")))
    tot_gene = tot_gene[!is.na(tot_gene)]

    N = length(tot_gene)
    df2 = tot_gene
    df1 = 1:N

    x <- as.integer(df1)
    y <- as.numeric(df2)

    p0 <- c(mean(y[which(x == 1)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_roary_hist <- function (sp) {
    input_path =file.path(mat_roary_hist_path, sp, "panmatrix.txt")
    panmatrix = read.csv(input_path, header=F, sep=" ")
    counts = table(rowSums(panmatrix))
    N = ncol(panmatrix)
    tab <- rep(0, N)
    names(tab) <- 1:N
    tab[names(counts)] <- counts
    tab
}
#--------------------------------------------------------------------------------
# PROKKA
#--------------------------------------------------------------------------------
read_prokka <- function (sp) {
    tot_gene = as.vector(unlist(read.csv(file.path(prokka_path, sp), header=F, sep=" ")))
    tot_gene = tot_gene[!is.na(tot_gene)]
    new_gene = vector("numeric", length(tot_gene))
    new_gene[1] = tot_gene[1]
    for (j in 2:(length(tot_gene))) {
        new_gene[j] = tot_gene[j] - tot_gene[j-1]
    }

    N = length(tot_gene)
    df2 = new_gene
    df1 = 1:N

    x <- as.integer(df1)
    y <- as.numeric(df2)
    # Remove x=1
    y <- y[x!=1]
    x <- x[x!=1]

    p0 <- c(mean(y[which(x == 2)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

#--------------------------------------------------------------------------------
get_k = function(len, a=4, p=0.0001) {
    return(ceiling(log((-len)/log(1-p),base=a)))
}
ceil_to_odd = function(num) {
    if (num%%2==0) {
        num = num + 1
    }
    return(num)
}
#--------------------------------------------------------------------------------
# KMER
#--------------------------------------------------------------------------------
read_kmer <- function (sp, k=NA) {
    if(is.na(k)) {
        k = ceil_to_odd(get_k(species_df() %>% filter(species==sp) %>% select(avg_sum_len)))
    }

    input_path =file.path(kmer_path, sp, paste0("pangenome_K:",k,"_can.txt"))
    tot_kmer = as.vector(unlist(read.csv(input_path, header=F, sep=" ")))
    tot_kmer = tot_kmer[!is.na(tot_kmer)]
    new_kmer = vector("numeric", length(tot_kmer))
    new_kmer[1] = tot_kmer[1]
    for (i in 2:(length(tot_kmer))) {
        new_kmer[i] = tot_kmer[i] - tot_kmer[i-1]
    }

    N = length(tot_kmer)
    df2 = new_kmer
    df1 = 1:N

    x <- as.integer(df1)
    y <- as.numeric(df2)
    # Remove x=1
    y <- y[x!=1]
    x <- x[x!=1]

    p0 <- c(mean(y[which(x == 2)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_kmer_tot <- function (sp, k=NA) {
    if(is.na(k)) {
        k = ceil_to_odd(get_k(species_df() %>% filter(species==sp) %>% select(avg_sum_len)))
    }
    input_path =file.path(kmer_path, sp, paste0("pangenome_K:",k,"_can.txt"))
    tot_kmer = as.vector(unlist(read.csv(input_path, header=F, sep=" ")))
    tot_kmer = tot_kmer[!is.na(tot_kmer)]

    N = length(tot_kmer)
    df2 = tot_kmer
    df1 = 1:N

    x <- as.integer(df1)
    y <- as.numeric(df2)

    p0 <- c(mean(y[which(x == 1)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_kmer_hist <- function (sp, k=NA) {
    if(is.na(k)) {
        k = ceil_to_odd(get_k(species_df() %>% filter(species==sp) %>% select(avg_sum_len)))
    }
    input_path =file.path(kmer_path, sp, paste0("hist_",k,"_can.txt"))
    h = as.vector(unlist(read.csv(input_path, header=F, sep=" ")))
    return(h)
}

#--------------------------------------------------------------------------------
# BPGA
#--------------------------------------------------------------------------------
# BPGA    
read_BPGA = function(filepath) {
    df = read.table(filepath, sep="\t", skip=1, fill=TRUE)
    df$V3 = NULL
    colnames(df) = c("x","y")
    df <- df[rowSums(is.na(df)) != ncol(df), ] # Remove row with all NAs
    return(df)
}

read_bpga_tot_mean <- function (sp, cluster="usearch", permutation=500, id=95) {
    input_path =file.path(bpga_path, cluster, paste0("perm",permutation,"_id",id), sp, "Supporting_files", "pan_genome.txt")
    df = read_BPGA(input_path)
    df = df %>% group_by(x) %>% summarise(y=mean(y)) # remove this to take all the value a not only the mean

    N = max(df$x)
    perm = sum(df$x==2)
    x <- as.integer(df$x)
    y <- as.numeric(df$y)

    p0 <- c(mean(y[which(x == 1)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_bpga_new <- function (sp, cluster="usearch", permutation=500, id=95, avg=TRUE) {
    input_path =file.path(bpga_path, cluster, paste0("perm",permutation,"_id",id), sp, "Supporting_files", "pan_genome.txt")
    df = read_BPGA(input_path)
    df = df %>% group_by(x) %>% summarise(y=mean(y)) 
    N = max(df$x)
    new_gene = vector("numeric", N)
    new_gene[1] = df$y[1]
    for (j in 2:N) {
        new_gene[j] = df$y[j] - df$y[j-1]
    }

    df2 = new_gene
    df1 = 1:N

    x <- as.integer(df1)
    y <- as.numeric(df2)
    # Remove x=1
    y <- y[x!=1]
    x <- x[x!=1]

    p0 <- c(mean(y[which(x == 2)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_mat_bpga_new<- function (sp, cluster="usearch", permutation=500, id=95) {
    cmd = paste("pangrowth growth -p",file.path(bpga_path, cluster, paste0("perm",permutation,"_id",id), 
                                                sp, "Supporting_files/matrix.txt"))
    growth = system(cmd, intern=T) 
    y = as.vector(sapply(strsplit(growth,split=" "),as.numeric))
    N = length(y)

    y_new = vector("numeric", N)
    y_new[1] = y[1]
    for (i in 2:(length(y))) {
        y_new[i] = y[i] - y[i-1]
    }

    x = 2:N
    y_new = y_new[x]

    p0 <- c(mean(y_new[which(x == 2)] ), 1)
    return(list(x=x,y=y_new,p0=p0))
}

read_mat_bpga_tot<- function (sp, cluster="usearch", permutation=500, id=95) {
    cmd = paste("pangrowth growth -p",file.path(bpga_path, cluster, paste0("perm",permutation,"_id",id), 
                                                sp, "Supporting_files/matrix.txt"))
    growth = system(cmd, intern=T) 
    y = as.vector(sapply(strsplit(growth,split=" "),as.numeric))
    N = length(y)
    x = 1:N

    p0 <- c(mean(y[which(x == 1)] ), 1)
    return(list(x=x,y=y,p0=p0))
}

read_bpga_hist <- function (sp, cluster="usearch", permutation=500, id=95) {
    input_path = file.path(bpga_path, cluster, paste0("perm",permutation,"_id",id), 
                           sp, "Supporting_files/matrix.txt")
    panmatrix = read.csv(input_path, header=F, sep="\t")
    counts = table(rowSums(panmatrix))
    N = ncol(panmatrix)
    tab <- rep(0, N)
    names(tab) <- 1:N
    tab[names(counts)] <- counts
    tab
}

#--------------------------------------------------------------------------------
# BENCHMARKS
#--------------------------------------------------------------------------------
get_benchmark <- function (program) {
    regex=NULL
    spec = c()
    benchmark_files = c()
    genome_num = c()
    tot_kmer = c()
    tot_bps = c()
    avg_sum_lens = c()
    genome_nums = c()

    for (sp in species()) {
        # Get the benchmark_files
        if (program=="kmer") {
            if(is.na(k)) {
                k = ceil_to_odd(get_k(species_df() %>% filter(species==sp) %>% select(avg_sum_len)))
            }
            regex = paste0("_k",k,".txt")
        }
        benchmarks_sp = list.files(file.path(folder_benchmarks, sp, program), pattern=regex, full.names = T)

        # Get info from seqkit
        seqkit = summarise_seqkit(file.path(folder_seqkit, sp, "fasta_stats.csv"))

        tot_bps = c(tot_bps, seqkit$tot_bp)
        avg_sum_lens = c(avg_sum_lens, seqkit$avg_sum_len)
        genome_nums = c(genome_nums, seqkit$N)

        spec = c(spec, sp)
        benchmark_files = c(benchmark_files, benchmarks_sp)
    }

    bfr = data.frame(matrix(nrow=0,ncol=0))
    for (bm in benchmark_files) {
        bfr = rbind(bfr, read.table(bm, header=T,sep="\t"))
    }

    bfr["species"] = spec
    bfr["tot_bp"] = tot_bps
    bfr["avg_sum_len"] = avg_sum_lens
    bfr["genome_num"] = genome_nums
    bfr
}

get_benchmark_prokka <- function () {
    spec = c()
    benchmark_files = c()

    for (sp in species()) {
        # Get the benchmark_files
        benchmarks_sp = list.files(file.path(folder_benchmarks, sp, "prokka"), full.names = T)
        spec = c(spec, rep(sp, length(benchmarks_sp)))
        benchmark_files = c(benchmark_files, benchmarks_sp)
    }

    bfr = data.frame(matrix(nrow=0,ncol=0))
    for (bm in benchmark_files) {
        bfr = rbind(bfr, read.table(bm, header=T,sep="\t"))
    }


    bfr["species"] = spec

    bfr = bfr %>% group_by(species) %>% 
        summarise(sum_s = sum(s), #time
                  max_rss=max(max_rss),max_uss=max(max_uss), max_pss=max(max_pss) #space
        )

    as.data.frame(bfr)
}

get_benchmark_bpga <- function (cluster="usearch", permutation=500, id=95) {
    spec = c()
    benchmark_files = c()
    genome_num = c()
    tot_kmer = c()
    tot_bps = c()
    avg_sum_lens = c()
    genome_nums = c()

    for (sp in species()) {
        # Get the benchmark_files
        benchmarks_sp = file.path(bpga_path, cluster, paste0("perm",permutation,"_id",id), sp, "time.txt")

        # Get info from seqkit
        seqkit = summarise_seqkit(file.path(folder_seqkit, sp, "fasta_stats.csv"))

        tot_bps = c(tot_bps, seqkit$tot_bp)
        avg_sum_lens = c(avg_sum_lens, seqkit$avg_sum_len)
        genome_nums = c(genome_nums, seqkit$N)

        spec = c(spec, sp)
        benchmark_files = c(benchmark_files, benchmarks_sp)
    }

    bfr = data.frame(matrix(nrow=0,ncol=0))
    for (bm in benchmark_files) {
        bfr = rbind(bfr, read_usrbintime(bm))
    }

    bfr["species"] = spec
    bfr["tot_bp"] = tot_bps
    bfr["avg_sum_len"] = avg_sum_lens
    bfr["genome_num"] = genome_nums
    bfr
}

# /usr/bin/time
# v stands for -v,verbos
read_usrbintime = function(filepath, type="v") {
    lines = readLines(filepath)

    kbytes_char = strsplit(lines[grep("Maximum resident set size",lines)],split=": ")[[1]][2]
    kbytes = as.numeric(kbytes_char)
    MRSS = kbytes/1024 #MB

    wall_clock_char = strsplit(lines[grep("Elapsed \\(wall clock\\) time",lines)],split=": ")[[1]][2]
    if(stringr::str_count(wall_clock_char, ":") == 1) {
        wall_clock = lubridate::period_to_seconds(lubridate::ms(wall_clock_char))
    } else if (stringr::str_count(wall_clock_char, ":") == 2) {
        wall_clock = lubridate::period_to_seconds(lubridate::hms(wall_clock_char))
    }
    data.frame(wall_clock=as.numeric(wall_clock), max_rss=as.numeric(MRSS))
}

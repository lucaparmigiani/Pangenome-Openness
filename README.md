# Pangenome Openness

This repository is dedicated to the experiment section of the paper:
Parmigiani, L., Wittler, R., Stoye, J.: [Revisiting pangenome openness with k-mers](https://doi.org/10.1101/2022.11.15.516472  ). bioRxiv. (2022).

- The analysis was performed on twelve bacterial pangenomes
- Tools:
  - k-mer-based approach
    - [Pangrowth](https://gitlab.ub.uni-bielefeld.de/gi/pangrowth )
  - gene-based approaches: 
    - [Roary](https://github.com/sanger-pathogens/Roary )
    - [Pantools](https://pantools.readthedocs.io/en/latest/index.html )  
    - [BPGA](https://iicb.res.in/bpga/ )
- Datasets: 12 bacterial species, from NCBI RefSeq, with the filter “Assembly level: Complete genome”. 
- Annotations: [Prokka](https://github.com/tseemann/prokka )(1.14.6, standard parameters). 

## Analysis
**Note**: some of these steps can take several hours (e.g., gene homology
clustering, annotation), therefore we have uploaded the final results of all the steps
with the scripts to work with them in the folder `analysis`. 

You can use directly that with the R scripts `analysis/reading.R` and
`analysis/fitting.R` to read and fit the data.

## Raw data
To generate the raw data from each we provide some scripts and a snakemake workflow.
The pipeline is divided into steps, so that it is not necessary to run it all if
you want to reproduce just some parts of the results. 

Each tool can be run using snakemake, except for BPGA for which a script is provided. 
The folder `config` contains snakemake config files that specify on which
bacteria the pipeline should be run and for which value of k. As before, three
options are given (`config_all`, `config_openclosed`, `config_quicktest`) but
any combination of species can be created, it is enough to provide the `kingdom`
(in this case `Bacteria` for all), `genus` (e.g., `Francisella`) and `species`
(e.g., `tularensis`).

## Table of Contents
<!--ts-->
   * [[Experiments] Revisiting pangenome openness with k-mers](#experiments-revisiting-pangenome-openness-with-k-mers)
      * [Analysis](#analysis)
      * [Raw data](#raw-data)
      * [Table of Contents](#table-of-contents)
      * [1. Data access](#1-data-access)
      * [2. Item extractions](#2-item-extractions)
         * [Pangrowth (k-mers)](#pangrowth-k-mers)
            * [Commands](#commands)
               * [kmer_run](#kmer_run)
         * [Genes](#genes)
            * [Annotation (Prokka)](#annotation-prokka)
            * [Roary](#roary)
               * [Dependencies](#dependencies)
               * [Commands](#commands-1)
                  * [roary_run](#roary_run)
            * [Pantools](#pantools)
               * [Dependencies](#dependencies-1)
               * [Commands](#commands-2)
                  * [pantools_run](#pantools_run)
            * [BPGA](#bpga)
      * [3. Analysis](#3-analysis)

<!-- Added by: luca, at: Thu Aug 17 01:42:08 PM CEST 2023 -->

<!--te-->

## 1. Data access

Genomes were downloaded from NCBI using [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download ).

```bash
pip install ncbi-genome-download
```

If the NCBI API changed or the tool is not working anymore, the folder
`data/accession` provides the accession numbers for each genome used.
These are also used to download the exact dataset we used.

In the folder we provide a script (`data/download.sh`) to download all the
genomes from NCBI. Note that some steps requires the genome to be uncompressed.
The script automatically uncompress the genomes and leaves them uncompressed for
the whole time. The whole uncompressed dataset is around 5GB.  
Each species will be downloaded in the folder `data`.

Usage:
```bash
./data/download.sh all
```

If you want to download the dataset used for the most of the images you can run:
```bash
./data/download.sh openclosed
```
which downloads only the genomes for _Streptococcus pneumoniae_ and _Yersinia pestis_.

For a quicker test there is also:
```bash
./data/download.sh quicktest
```
Which downloads only 57 genomes of _Francisella tularensis_.

## 2. Item extractions
### Pangrowth (k-mers)

- [Pangrowth](https://gitlab.ub.uni-bielefeld.de/gi/pangrowth )

**Install pangrowth** 
```bash
cd scripts
git clone https://gitlab.ub.uni-bielefeld.de/gi/pangrowth
cd pangrowth 
make
cd ../..
```
#### Commands

##### kmer_run

This produce the histogram and the pangenome growth for **k-mers** for each species
in the `config` file. The output can be find in `results/data/kmer_k/species` 
(e.g., `results/data/kmer_17/Francisella_tularensis`).

```bash
snakemake --cores 12 --latency-wait 60 kmer_run -p --verbose --rerun-incomplete --configfile config/config_quicktest.yaml
snakemake --cores 12 --latency-wait 60 kmer_run -p --verbose --rerun-incomplete --configfile config/config_openclosed.yaml
snakemake --cores 12 --latency-wait 60 kmer_run -p --verbose --rerun-incomplete --configfile config/config_all.yaml
```

### Genes

#### Annotation (Prokka)

- [Prokka](https://github.com/tseemann/prokka )

The genes pipeline requires the annotation with Prokka. This step can be
run manually or it will be called automatically by running either `roary_run`
or `pantools_run`. The Annotations are also required if you want to run BPGA.

#### Roary

##### Dependencies

- [Roary](https://github.com/sanger-pathogens/Roary )

**Ubuntu**
```bash
sudo apt-get install roary
```
- [Prokka](https://github.com/tseemann/prokka )

##### Commands

###### roary_run

This produce the histogram and the pangenome growth for **genes** for each species
in the `config` file. The output can be find in `results/data/gene_roary/species` 
(e.g., `results/data/gene_roary/Francisella_tularensis`).

```bash
snakemake --cores 12 --latency-wait 60 roary_run -p --verbose --rerun-incomplete --configfile config/config_quicktest.yaml
snakemake --cores 12 --latency-wait 60 roary_run -p --verbose --rerun-incomplete --configfile config/config_openclosed.yaml
snakemake --cores 12 --latency-wait 60 roary_run -p --verbose --rerun-incomplete --configfile config/config_all.yaml
```

#### Pantools

##### Dependencies

- [Pantools](https://pantools.readthedocs.io/en/latest/index.html )  
  * Refer to the
    [manual](https://pantools.readthedocs.io/en/latest/user_guide/install.html )
    for their intallations. 
  * Pantools must be installed in the directory `scripts/pantools`

**Ubuntu**
```bash
mkdir -p scripts/pantools
cd scripts/pantools
wget https://www.bioinformatics.nl/pangenomics/data/pantools-4.2.2.jar  
```
- [Prokka](https://github.com/tseemann/prokka )

##### Commands
###### pantools_run

This produce the histogram and the pangenome growth for **genes** for each species
in the `config` file. The output can be find in `results/data/gene_pantools/species` 
(e.g., `results/data/gene_pantools/Francisella_tularensis`).

```bash
snakemake --cores 12 --latency-wait 60 pantools_run -p --verbose --rerun-incomplete --configfile config/config_quicktest.yaml
snakemake --cores 12 --latency-wait 60 pantools_run -p --verbose --rerun-incomplete --configfile config/config_openclosed.yaml
snakemake --cores 12 --latency-wait 60 pantools_run -p --verbose --rerun-incomplete --configfile config/config_all.yaml
```

#### BPGA

BPGA can be installed from [https://sourceforge.net/projects/bpgatool/files/](https://sourceforge.net/projects/bpgatool/files/).
You can find a script at `scripts/bpga_automatic.sh` that accepts a species as
input and initiates a tmux session, automatically sending the necessary
keystrokes to run BPGA. 
Please ensure the `faa` variable within `scripts/bpga_automatic.sh` is set to
the directory containing all the relevant protein files.  
The script has to be copied into a directory containg BPGA as follows:

```bash
wget https://downloads.sourceforge.net/project/bpgatool/BPGA-1.3-linux-x86_64-0-0-0.tar.gz
tar -xf BPGA-1.3-linux-x86_64-0-0-0/
cp ./scripts/bpga_automatic.sh BPGA-1.3-linux-x86_64-0-0-0/BPGA-Version-1.3/bin
cd BPGA-1.3-linux-x86_64-0-0-0/BPGA-Version-1.3/bin
chmod +x BPGA-Version-1.3
```

If you encounter issues with the automated script, BPGA can also be executed
manually. Simply pass in the directory containing the protein files.

The folder containing all the proteins can be obtained by running Prokka.

## 3. Analysis

The folder `analysis` contains all the scripts used to compare the results of
the three methods and the resulting data. 
Since the computation of all the steps before can be onerous
we provided intermediate steps in the folder `analysis/data`.

- The file `analysis/pankmer_reading.R` provides functions to read the results from
  Pangrowth, Roary, Pantools, BPGA and Prokka.

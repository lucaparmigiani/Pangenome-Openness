# folder, name, ext, file, path
def ls(folder, regex=r'.', with_folder=False, with_ext=True):
    import re
    fastqre = re.compile(regex)

    if os.path.isdir(folder):
        listfile = list( map( lambda x: os.path.join(folder, x),  \
                          filter(fastqre.search,  \
                                 os.listdir(os.path.join(folder) ) ) ) )
        if listfile:
            if not with_folder:
                listfile = [os.path.basename(i) for i in listfile]
            if not with_ext:
                listfile = [os.path.splitext(i)[0] for i in listfile]

        return listfile
    else:
        print("folder:", folder,"does not exists")
        return []

#--------------------------------------------------------------------------------
def fasta_name(wildcards):
    return ls(os.path.join(DATA, wildcards.species), r'\.fna', False, False)

def fasta_path(wildcards):
    return expand(DATA+wildcards.species+"/{sample}.fna", sample=fasta_name(wildcards))
    
def gff_path(wildcards):
    return expand(OUTDIR+"prokka/"+wildcards.species+"/{sample}.gff", sample=fasta_name(wildcards))

def get_genus(wildcards):
    return config["bacteria"][wildcards.species]["genus"]

def get_kingdom(wildcards):
    return config["bacteria"][wildcards.species]["kingdom"]

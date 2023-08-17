#pip install --user ncbi-genome-download

source=${BASH_SOURCE[0]}

while [ -L "$source" ]; do 
    owndir=$( cd -P "$( dirname "$source" )" >/dev/null 2>&1 && pwd )
    source=$(readlink "$source")
    [[ $source != /* ]] && source=$owndir/$source 
done

owndir=$( cd -P "$( dirname "$source" )" >/dev/null 2>&1 && pwd )

accession_dir=$owndir/accession
outdir=$owndir

function show_usage() {
    echo usage:
    echo ./download.sh all
    echo "(download all the 12 species.)"
    echo
    echo ./download.sh openclosed
    echo "(download only Streptococcus pneumoniae and Yersinia pestis.)"
    echo
    echo ./download.sh quicktest
    echo "(download only 8 genomes of Francisella tularensis.)"
}

function unzip() {
    species=$1
    echo Uncompressing $outdir/$species
    mv $(du -aL $outdir/$species | grep gz$ | cut -f2) $outdir/$species
    gzip -df  $outdir/$species/*.gz
    rm -r $outdir/$species/refseq
}

if [[ "$#" == 1 ]]; then
    if [[ $1 == "openclosed" ]]; then
        echo Downloading $outdir/Yersinia_pestis
        ncbi-genome-download -s refseq --assembly-levels complete --formats fasta -A ${accession_dir}/Yersinia_pestis.txt bacteria -p 4 -r -4 -o $outdir/Yersinia_pestis
        unzip Yersinia_pestis
        echo Downloading $outdir/Streptococcus_pneumoniae
        ncbi-genome-download -s refseq --assembly-levels complete --formats fasta -A ${accession_dir}/Streptococcus_pneumoniae.txt bacteria -p 4 -r -4 -o $outdir/Streptococcus_pneumoniae
        unzip Streptococcus_pneumoniae
        exit

    elif [[ $1 == "all" ]]; then
        for acc in $(ls ${accession_dir}/*.txt); do
            species=$(basename -- "$acc")
            species="${species%.*}"
            echo Downloading $outdir/$species
            ncbi-genome-download -s refseq --assembly-levels complete --formats fasta -A $acc bacteria -p 4 -r -4 -o $outdir/$species
            unzip $species
        done
        exit

    elif [[ $1 == "quicktest" ]]; then
        echo Downloading $outdir/Francisella_tularensis
        ncbi-genome-download -s refseq --assembly-levels complete --formats fasta -A ${accession_dir}/Francisella_tularensis.txt bacteria -p 4 -r -4 -o $outdir/Francisella_tularensis
        unzip Francisella_tularensis
        exit

    else
        echo $1 not recognized
    fi

else
    echo Error: only one parameter accepted
fi

show_usage

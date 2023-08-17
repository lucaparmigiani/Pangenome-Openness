#possible species: 
#"Bacillus_cereus" "Buchnera_aphidicola" "Campylobacter_jejuni"
#"Clostridium_botulinum" "Coxiella_burnetii" "Francisella_tularensis"
#"Helicobacter_pylori" "Prochlorococcus_marinus" "Rhodopseudomonas_palustris"
#"Streptococcus_pneumoniae" "Streptococcus_pyogenes" "Yersinia_pestis"

species=$1
mkdir -p ../$species

faa=/vol/data/PanKmer/results/prokka/$species/faa
server=Server$species

tmux new-session -d -s $server 'bash'
#tmux send-keys -t $server -l "/usr/bin/time -o ../time_a.txt -f 'mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U' ./BPGA-Version-1.3"
tmux send-keys -t $server -l "/usr/bin/time -o ../$species/time.txt -v ./BPGA-Version-1.3"
tmux send-keys -t $server Enter
#tmux send-keys -t $server -l "./BPGA-Version-1.3"
#tmux send-keys -t $server Enter
tmux send-keys -t $server -l 1 #preparer fasta files
tmux send-keys -t $server Enter
tmux send-keys -t $server -l 4 #normal protein fasta
tmux send-keys -t $server Enter
tmux send-keys -t $server -l $faa
tmux send-keys -t $server Enter
tmux send-keys -t $server Enter
tmux send-keys -t $server -l 2 #start analysis
tmux send-keys -t $server Enter
tmux send-keys -t $server -l 1 #usearch
tmux send-keys -t $server Enter
tmux send-keys -t $server -l 0.95 #cutoff seq identity
tmux send-keys -t $server Enter
tmux send-keys -t $server -l 1 #choose for the pan-core drawing needed for the permutation step
tmux send-keys -t $server Enter
tmux send-keys -t $server -l 500 #permutations
tmux send-keys -t $server Enter
tmux send-keys -t $server -l 0 #exit
tmux send-keys -t $server Enter

tmux send-keys -t $server "cp -fr Results ../$species"
tmux send-keys -t $server Enter
tmux send-keys -t $server "cp -fr Supporting_files ../$species"
tmux send-keys -t $server Enter
tmux send-keys -t $server -l exit
tmux send-keys -t $server Enter

mkdir -p growth
for sp in $(ls -l | grep ^d | awk '{print $9}'); do
    echo $sp
    pangrowth growth -p $sp/panmatrix.txt > growth/$sp
done


for sp in $(ls | grep \.txt); do
    echo $sp
    /home/luca/@focus/soft/pangrowth/pangrowth growth -p -i $sp > growth/$sp
done


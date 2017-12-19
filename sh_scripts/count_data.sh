for i in $(find AVHRR*/ -maxdepth 0 -type d) ; do 
    echo -n $i": " ; 
    (find $i -type f -name '*.data' | wc -l) ; 
done


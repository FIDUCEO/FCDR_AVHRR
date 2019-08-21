for i in $(find AVHRR*/ -maxdepth 2 -type d) ; do 
    echo -n $i": " ; 
    (find $i -type f -name '*.nc' | wc -l) ; 
done



#!/usr/bin/env bash
# Run performance tests to compare the different implementations

NODES=4
#N="64 128 256 512"
#N=$(seq 64 4 512)
#N=$(nawk 'BEGIN{ for(i=64;i<=512;i=i+4) print i}')
N=$(nawk 'BEGIN{ for(i=64;i<=1024;i=i+1) print i}')
OUTPUTFILE="tests_perf.tsv tests_perf2.tsv tests_perf3.tsv tests_perf4.tsv"

for OUTFILE in $OUTPUTFILE
do
    # Sanity checks
    [ -f $OUTFILE ] && echo "$OUTFILE file already exists, not overwriting (to be safe)!" && exit 1;

    echo "Performance tests:"
    for n in $N
    do
        echo "Running for $n"
        #S1=$(./yee -n $n                        | grep 'Elapsed' | awk '{print $2}')
        S2=$(./yee_omp -n $n -t $NODES          | grep 'Elapsed' | awk '{print $2}')
        S3=$(./yee_pthr -n $n -t $NODES         | grep 'Elapsed' | awk '{print $2}')
        S4=$(mpirun -np $NODES ./yee_mpi -n $n  | grep 'Elapsed' | awk '{print $2}')
        S5=$(mpirun -np $NODES ./yee_mpi2 -n $n | grep 'Elapsed' | awk '{print $2}')
        echo $n $S2 $S3 $S4 $S5 >> $OUTFILE
    done
done


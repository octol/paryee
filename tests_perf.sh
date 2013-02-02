#!/usr/bin/env bash
# Run performance tests to compare the different implementations

NNODES="4 8"
#N=$(seq 64 4 512)
N=$(nawk 'BEGIN{ for(i=100;i<=1000;i=i+1) print i}')
OUTPUTFILE="tests_perf1 tests_perf2 tests_perf3 tests_perf4"

for NODES in $NNODES
do
    for OUTFILE in $OUTPUTFILE
    do
        OUTFILE=${OUTFILE}_${NODES}.tsv

        # Sanity checks
        [ -f $OUTFILE ] && echo "$OUTFILE file already exists, not overwriting (to be safe)!" && exit 1;

        echo "Performance tests (${NODES} nodes):"
        for n in $N
        do
            echo "Running for $n"
            S1=$(./yee_omp -q -n $n -t $NODES          | grep 'Elapsed' | awk '{print $2}')
            S2=$(./yee_pthr -q -n $n -t $NODES         | grep 'Elapsed' | awk '{print $2}')
            S3=$(mpirun -np $NODES ./yee_mpi -q -n $n  | grep 'Elapsed' | awk '{print $2}')
            S4=$(mpirun -np $NODES ./yee_mpi2 -q -n $n | grep 'Elapsed' | awk '{print $2}')
            echo $n $S1 $S2 $S3 $S4 >> $OUTFILE
        done
    done
done



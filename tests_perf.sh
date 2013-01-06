#!/usr/bin/env bash
# Run performance tests to compare the different implementations

NODES=4
#N="64 128 256 512"
#N=$(seq 64 4 512)
N=$(nawk 'BEGIN{ for(i=64;i<=512;i=i+4) print i}')
OUTFILE=tests_perf.tsv

# Sanity checks
[ -f $OUTFILE ] && echo "$OUTFILE file already exists, not overwriting (to be safe)!" && exit 1;

echo "Performance tests:"
for n in $N
do
    echo "Running for $n"
    S1=$(./yee -n $n                       | grep 'Elapsed' | awk '{print $2}')
    S2=$(./yee_omp -n $n -t $NODES         | grep 'Elapsed' | awk '{print $2}')
    S3=$(./yee_pthr -n $n -t $NODES        | grep 'Elapsed' | awk '{print $2}')
    S4=$(mpirun -np $NODES ./yee_mpi -n $n | grep 'Elapsed' | awk '{print $2}')
    echo $n $S1 $S2 $S3 $S4 >> $OUTFILE
done


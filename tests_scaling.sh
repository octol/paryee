#!/usr/bin/env bash
# Run scaling tests to compare the different implementations

#NODES=$(seq 8)
#NODES="1 2 4 8"
NODES=$(nawk 'BEGIN{ for(i=1;i<=8;i++) print i}')
N=512
OUTFILE=tests_scaling.tsv

# Sanity checks
[ -f $OUTFILE ] && echo "$OUTFILE file already exists, not overwriting (to be safe)!" && exit 1;

echo "Scaling tests:"
for M in $NODES
do
    echo "Running for $M threads"
    S1=$(./yee -n $N                   | grep 'Elapsed' | awk '{print $2}')
    S2=$(./yee_omp -n $N -t $M         | grep 'Elapsed' | awk '{print $2}')
    S3=$(./yee_pthr -n $N -t $M        | grep 'Elapsed' | awk '{print $2}')
    S4=$(mpirun -np $M ./yee_mpi -n $N | grep 'Elapsed' | awk '{print $2}')
    echo $M $S1 $S2 $S3 $S4 >> $OUTFILE
done

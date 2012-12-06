#!/usr/bin/env bash
# Run scaling tests to compare the different implementations

NODES="1 2 4 8"
#N=1024
#N=8192
N=32768
#N=65536
OUTFILE=tests_scaling.tsv

# Sanity checks
[ -f $OUTFILE ] && echo "$OUTFILE file already exists, not overwriting (to be safe)!" && exit 1;

echo "Scaling tests:"
for M in $NODES
do
    MPI_NODES=$((M + 1))
    echo "Running for $M threads"
    S1=$(./yee -n $N                          | grep 'Elapsed' | awk '{print $2}')
    S2=$(./yee_ref -n $N                      | grep 'Elapsed' | awk '{print $2}')
    S3=$(./yee_omp -n $N -t $M                | grep 'Elapsed' | awk '{print $2}')
    S4=$(./yee_pthr -n $N -t $M               | grep 'Elapsed' | awk '{print $2}')
    S5=$(mpirun -n $MPI_NODES ./yee_mpi -n $N | grep 'Elapsed' | awk '{print $2}')
    echo $M $S1 $S2 $S3 $S4 $S5 >> $OUTFILE
done

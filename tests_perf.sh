#!/usr/bin/env bash
# Run performance tests to compare the different implementations

NODES=4
NN="1024 2048 4096 8192 16384"
#NN="1024 2048 4096 8192 16384 32768"
#NN="8192 16384 32768 65536"
OUTFILE=tests_perf.tsv

MPI_NODES=$((NODES + 1))

# Sanity checks
[ -f $OUTFILE ] && echo "$OUTFILE file already exists, not overwriting (to be safe)!" && exit 1;

echo "Performance tests:"
for N in $NN
do
    echo "Running for $N"
    S1=$(./yee -n $N                          | grep 'Elapsed' | awk '{print $2}')
    S2=$(./yee_ref -n $N                      | grep 'Elapsed' | awk '{print $2}')
    S3=$(./yee_omp -n $N -t $NODES            | grep 'Elapsed' | awk '{print $2}')
    S4=$(./yee_pthr -n $N -t $NODES           | grep 'Elapsed' | awk '{print $2}')
    S5=$(mpirun -n $MPI_NODES ./yee_mpi -n $N | grep 'Elapsed' | awk '{print $2}')
    echo $N $S1 $S2 $S3 $S4 $S5 >> $OUTFILE
done


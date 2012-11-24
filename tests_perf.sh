#!/usr/bin/env bash
# Run performance tests to compare the different implementations

NODES=2
NN="256 512 1024 2048 4096 8162"
OUTFILE=tests_perf.tsv

MPI_NODES=$((NODES + 1))

# Sanity checks
if [ -f $OUTFILE ]
then
    echo "Output file already exists!"
    exit 1;
fi

# Run tests
for N in $NN
do
    echo "Running for $N"
    S1=$(./yee -n $N                         | grep 'Elapsed' | awk '{print $2}')
    S2=$(./yee_ref -n $N                     | grep 'Elapsed' | awk '{print $2}')
    S3=$(./yee_omp -n $N -t $NODES           | grep 'Elapsed' | awk '{print $2}')
    S4=$(./yee_pthr -n $N -t $NODES          | grep 'Elapsed' | awk '{print $2}')
    S5=$(./yee_pthr_barrier -n $N -t $NODES  | grep 'Elapsed' | awk '{print $2}')
    S6=$(./yee_mpi -n $N -t $MPI_NODES       | grep 'Elapsed' | awk '{print $2}')
    echo $N $S1 $S2 $S3 $S4 $S5 $S6 >> $OUTFILE
done


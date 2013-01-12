#!/usr/bin/env bash
# Run scaling tests to compare the different implementations

#NODES=$(seq 8)
NODES=$(nawk 'BEGIN{ for(i=1;i<=8;i++) print i}')
N=500
OUTPUTFILE="tests_scaling1 tests_scaling2 tests_scaling3 tests_scaling4"

for OUTFILE in $OUTPUTFILE
do
    OUTFILE=${OUTFILE}.tsv

    # Sanity checks
    [ -f $OUTFILE ] && echo "$OUTFILE file already exists, not overwriting (to be safe)!" && exit 1;

    echo "Scaling tests: ${OUTFILE}"
    for M in $NODES
    do
        echo "Running for $M threads"
        #S1=$(./yee -n $N                    | grep 'Elapsed' | awk '{print $2}')
        S2=$(./yee_omp -n $N -t $M          | grep 'Elapsed' | awk '{print $2}')
        S3=$(./yee_pthr -n $N -t $M         | grep 'Elapsed' | awk '{print $2}')
        S4=$(mpirun -np $M ./yee_mpi -n $N  | grep 'Elapsed' | awk '{print $2}')
        S5=$(mpirun -np $M ./yee_mpi2 -n $N | grep 'Elapsed' | awk '{print $2}')
        echo $M $S2 $S3 $S4 $S5 >> $OUTFILE
    done
done

# Now merge computed data and write to tests_perf.tsv
OUTFILE_M=tests_scaling.tsv
[ -f $OUTFILE_M ] && echo "$OUTFILE_M file already exists, not overwriting (to be safe)!" && exit 1;
octave -qf gather_scaling_data.m > ${OUTFILE_M}


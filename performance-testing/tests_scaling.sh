#!/usr/bin/env bash
# Run scaling tests to compare the different implementations

bindir=bin
outdir=output-tests

#nodes=$(seq 8)
nodes=$(nawk 'BEGIN{ for(i=1;i<=8;i++) print i}')
#N=1000
N=32
outputfile="tests_scaling1 tests_scaling2 tests_scaling3 tests_scaling4"

if [ ! -f README.md ]; then
    echo "Please run from project root directory"
    exit 1;
fi

for outfile in ${outputfile}
do
    outfile=${outdir}/${outfile}.tsv

    # Sanity checks
    if [ -f ${outfile} ]; then
        echo "${outfile} file already exists, not overwriting!" 
        exit 1;
    fi

    echo "Scaling tests: ${outfile}"
    for M in ${nodes}
    do
        echo "Running for $M threads"
        #S1=$(./${bindir}/yee -n $N                       | grep 'Elapsed' | awk '{print $2}')
        S2=$(./${bindir}/yee_omp -n $N -t $M -q          | grep 'Elapsed' | awk '{print $2}')
        S3=$(./${bindir}/yee_pthr -n $N -t $M -q         | grep 'Elapsed' | awk '{print $2}')
        S4=$(mpirun -np $M ./${bindir}/yee_mpi -n $N -q  | grep 'Elapsed' | awk '{print $2}')
        S5=$(mpirun -np $M ./${bindir}/yee_mpi2 -n $N -q | grep 'Elapsed' | awk '{print $2}')
        echo $M $S2 $S3 $S4 $S5 >> ${outfile}
    done
done


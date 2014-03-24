#!/usr/bin/env bash
# Run performance tests to compare the different implementations

bindir=bin
outdir=output-tests

nnodes="4 8"
#N=$(seq 64 4 512)
N=$(nawk 'BEGIN{ for(i=100;i<=1000;i=i+1) print i}')
outputfile="tests_perf1 tests_perf2 tests_perf3 tests_perf4"

if [ ! -f README.md ]; then
    echo "Please run from project root directory"
    exit 1;
fi

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

for nodes in ${nnodes}; do
    for outfile in ${outputfile}; do
        outfile=${outdir}/${outfile}_${nodes}.tsv

        # Sanity checks
        if [ -f ${outfile} ]; then
            echo "${outfile} file already exists, not overwriting!" 
            exit 1;
        fi

        echo "Performance tests (${nodes} nodes):"
        for n in $N; do
            echo "Running for $n"
            S1=$(./${bindir}/yee_omp -q -n $n -t ${nodes}          | grep 'Elapsed' | awk '{print $2}')
            S2=$(./${bindir}/yee_pthr -q -n $n -t ${nodes}         | grep 'Elapsed' | awk '{print $2}')
            S3=$(mpirun -np ${nodes} ./${bindir}/yee_mpi -q -n $n  | grep 'Elapsed' | awk '{print $2}')
            S4=$(mpirun -np ${nodes} ./${bindir}/yee_mpi2 -q -n $n | grep 'Elapsed' | awk '{print $2}')
            echo $n $S1 $S2 $S3 $S4 >> ${outfile}
        done
    done
done



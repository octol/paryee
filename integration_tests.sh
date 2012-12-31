#!/usr/bin/env bash
# Run integration tests by comparing to single thread reference solution.

N=64
threads=4
limit=1e-14

yee_bin="yee_omp yee_pthr yee_mpi"
args="-n $N" 

printf "\n  Integration tests:\n"
printf "  Testing that all implementations produce the same output.\n\n"

# Generate reference solution to compare against
./yee -o /tmp/yee.tsv $args 1> /dev/null

for yb in $yee_bin
do
    # Compute solution
    outfile=/tmp/${yb}.tsv
    if [ $yb == "yee_mpi" ]
    then
        mpirun -n $threads ./$yb -o $outfile $args 1> /dev/null
    else
        ./$yb -t $threads -o $outfile $args 1> /dev/null
    fi

    # Compare
    errlog=/tmp/numdiff_log_${yb}.txt
    numdiff -S -a $limit /tmp/yee.tsv $outfile > $errlog

    # Eval results
    if [ $? == 0 ]
    then
        echo "OK: $yb passed test."
    else
        echo "FAIL: $yb output differs! See $errlog"
    fi
done



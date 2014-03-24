#!/usr/bin/env bash
# Run integration tests by comparing to single thread reference solution.

N=64
threads=4
limit=1e-14

yee_bin="yee_omp yee_pthr yee_mpi yee_mpi2"
args="-n $N" 

bindir=bin

printf "\n  Integration tests:\n"
printf "  Testing that all implementations produce the same output.\n"

# Generate reference solution to compare against
./$bindir/yee -o /tmp/yee.tsv $args 1> /dev/null

for yb in $yee_bin
do
    printf "\nTest: ${yb}\n"
    # Compute solution
    outfile=/tmp/${yb}.tsv
    if [ $yb == "yee_mpi" ] || [ $yb == "yee_mpi2" ]
    then
        mpirun -n $threads ./$bindir/$yb -o $outfile $args 1> /dev/null
    else
        ./$bindir/$yb -t $threads -o $outfile $args 1> /dev/null
    fi

    # Compare
    errlog=/tmp/numdiff_log_${yb}.txt
    numdiff -S -a $limit /tmp/yee.tsv $outfile > $errlog
    #numdiff -S -a $limit /tmp/yee.tsv $outfile | grep -A 1 "Largest absolute error"
    success=$?
    cat $errlog | grep -A 1 "Sum of all absolute errors"

    # Eval results
    if [ $success == 0 ]
    then
        echo "OK: $yb passed test."
        diff /tmp/yee.tsv $outfile > $errlog
        if [ $? == 0 ]
        then
            echo "Results exactly identical."
        fi
    else
        echo "FAIL: $yb output differs! See $errlog"
    fi
done



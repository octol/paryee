#!/usr/bin/env bash
# Run integrity checks by comparing to single thread reference solution.

N=128
M=4
M_MPI=5

yee_bin="yee yee_omp yee_pthr yee_mpi"

# Generate reference solution to compare against
./yee_ref -n $N -p /tmp/yee_ref_p.tsv -u /tmp/yee_ref_u.tsv

for yb in $yee_bin
do
    # Compute solution
    echo "*** Testing: $yb"
    if [ $yb == "yee_mpi" ]
    then
        mpirun -n $M_MPI ./$yb -n $N -p /tmp/yee_p.tsv -u /tmp/yee_u.tsv
    else
        ./$yb -n $N -p /tmp/yee_p.tsv -u /tmp/yee_u.tsv
    fi

    # Compare
    ERRORLOG=/tmp/numdiff_log_${yb}.txt
    numdiff -S -a 1e-6 /tmp/yee_ref_p.tsv /tmp/yee_p.tsv > $ERRORLOG

    # Eval results
    if [ $? == 0 ]
    then
        echo "OK: $yb passed test."
    else
        echo "FAIL: $yb output differs from reference implementation. See $ERRORLOG"
    fi
done



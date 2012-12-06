#!/usr/bin/env bash
# Run integrity checks by comparing to single thread reference solution.

N=1024
M=4
M_MPI=5

yee_bin="yee yee_omp yee_pthr yee_mpi"

# Generate reference solution to compare against
./yee_ref -n $N -p /tmp/yee_ref_p.tsv -u /tmp/yee_ref_u.tsv

for yb in $yee_bin
do
    # Compute solution
    echo " --- Testing: $yb"
    outfile_u=/tmp/${yb}_u.tsv
    outfile_p=/tmp/${yb}_p.tsv
    if [ $yb == "yee" ]
    then
        ./$yb -n $N -p $outfile_p -u $outfile_u
    elif [ $yb == "yee_ref" ]
    then
        ./$yb -n $N -p $outfile_p -u $outfile_u
    elif [ $yb == "yee_mpi" ]
    then
        mpirun -n $M_MPI ./$yb -n $N -p $outfile_p -u $outfile_u
    else
        ./$yb -n $N -t $M -p $outfile_p -u $outfile_u
    fi

    # Compare
    errlog=/tmp/numdiff_log_${yb}.txt
    numdiff -S -a 1e-7 /tmp/yee_ref_u.tsv $outfile_u > $errlog

    # Eval results
    if [ $? == 0 ]
    then
        echo "OK: $yb passed test."
    else
        echo "FAIL: $yb output differs from ref implementation. See $errlog"
    fi
done



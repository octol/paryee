#!/usr/bin/env bash
# Run performance tests to compare the different implementations

readonly PROGNAME=$(basename $0)
readonly PROGDIR=$(readlink -m $(dirname $0))
readonly ARGS="$@"

readonly HOSTNAME=$(hostname)
readonly ARCHNAME=$(uname -s)_$(uname -m)
readonly BINDIR="$ARCHNAME/bin"
readonly OUTDIR="perf-test-out/$HOSTNAME"

# Default arguments
NODES=4
N_MAX=1000
SAMPLES=4
OUTPUT_BASENAME="test_perf"
FORCE=no

cmdline() {
    while getopts "n:N:s:f" OPTION; do
        case "$OPTION" in
            n)  NODES="$OPTARG" ;;
            N)  N_MAX="$OPTARG" ;;
            s)  SAMPLES="$OPTARG" ;;
            f)  FORCE=yes ;;
            *)  echo Usage: "$PROGNAME": [-n NODES] [-N MAX_GRIDLENGTH] [-s SAMPLES] [-f]
                exit 1;;
        esac
    done
}

verify_environment() {
    [ -f README.md ] || ( echo "Run from project root dir"; exit 1 )
    [ -d $BINDIR ] || ( echo "$BINDIR: not found"; exit 1 )
    mkdir -p $OUTDIR
}

run_test() {
    local grid_length=$1
    local output_file=$2
    S1=$(./${BINDIR}/yee_omp -q -n $grid_length -t ${NODES}          | grep 'Elapsed' | awk '{print $2}')
    S2=$(./${BINDIR}/yee_pthr -q -n $grid_length -t ${NODES}         | grep 'Elapsed' | awk '{print $2}')
    S3=$(mpirun -np ${NODES} ./${BINDIR}/yee_mpi -q -n $grid_length  | grep 'Elapsed' | awk '{print $2}')
    S4=$(mpirun -np ${NODES} ./${BINDIR}/yee_nonblock_mpi -q -n $grid_length | grep 'Elapsed' | awk '{print $2}')
    echo $grid_length $S1 $S2 $S3 $S4 >> ${output_file}
}

run_all_tests() {
    local N=$(for ((i=100;i<=$N_MAX;i++)); do echo $i; done)

    for ((i=1; i <= $SAMPLES; i++)); do
        local output_file="${OUTDIR}/${OUTPUT_BASENAME}_${NODES}_${i}.tsv"

        if [[ -f ${output_file} && ! $FORCE == "yes" ]]; then
            echo "${output_file} exists, skipping!"
        else
            echo "Performance tests ($NODES nodes):"
            for grid_length in $N; do
                echo "Running for $grid_length"
                run_test $grid_length $output_file
            done
        fi
    done
}

main() {
    cmdline $ARGS
    verify_environment
    run_all_tests
}
main $ARGS

# nnodes="4 8"
# #N=$(seq 64 4 512)
# N=$(nawk 'BEGIN{ for(i=100;i<=1000;i=i+1) print i}')
# outputfile="tests_perf1 tests_perf2 tests_perf3 tests_perf4"
# 
# for nodes in ${nnodes}; do
#     for outfile in ${outputfile}; do
#         outfile=${outdir}/${outfile}_${nodes}.tsv
# 
#         echo "Performance tests (${nodes} nodes):"
#         for n in $N; do
#             echo "Running for $n"
#             S1=$(./${bindir}/yee_omp -q -n $n -t ${nodes}          | grep 'Elapsed' | awk '{print $2}')
#             S2=$(./${bindir}/yee_pthr -q -n $n -t ${nodes}         | grep 'Elapsed' | awk '{print $2}')
#             S3=$(mpirun -np ${nodes} ./${bindir}/yee_mpi -q -n $n  | grep 'Elapsed' | awk '{print $2}')
#             S4=$(mpirun -np ${nodes} ./${bindir}/yee_mpi2 -q -n $n | grep 'Elapsed' | awk '{print $2}')
#             echo $n $S1 $S2 $S3 $S4 >> ${outfile}
#         done
#     done
# done
# 
# 

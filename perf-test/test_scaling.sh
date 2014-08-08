#!/usr/bin/env bash
# Run scaling tests to compare the different implementations
#
# Example: Run 8 nodes on a 128x128 grid and sample 4 times
#
#   ./test_scaling.sh -n 8 -N 128 -s 4
#

readonly PROGNAME=$(basename $0)
readonly PROGDIR=$(readlink -m $(dirname $0))
readonly ARGS="$@"

readonly HOSTNAME=$(hostname)
readonly ARCHNAME=$(uname -s)_$(uname -m)
readonly BINDIR="$ARCHNAME/bin"
readonly OUTDIR="perf-test-out/$HOSTNAME"

# Default arguments
NUM_NODES=8
N=1000
SAMPLES=4
OUTPUT_BASENAME="test_scaling"
FORCE=no

cmdline() {
    while getopts "n:N:s:f" OPTION; do
        case "$OPTION" in
            n)  NUM_NODES="$OPTARG" ;;
            N)  N="$OPTARG" ;;
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
    local nodes="$1"
    local output_file="$2"
    S1=$(./${BINDIR}/yee_omp -n $N -t $nodes -q          | grep 'Elapsed' | awk '{print $2}')
    S2=$(./${BINDIR}/yee_pthr -n $N -t $nodes -q         | grep 'Elapsed' | awk '{print $2}')
    S3=$(mpirun -np $nodes ./${BINDIR}/yee_mpi -n $N -q  | grep 'Elapsed' | awk '{print $2}')
    S4=$(mpirun -np $nodes ./${BINDIR}/yee_nonblock_mpi -n $N -q | grep 'Elapsed' | awk '{print $2}')
    echo $nodes $S1 $S2 $S3 $S4 >> ${output_file}
}

run_all_tests() {
    local nodes=$(for ((i=1;i<=$NUM_NODES;i++)); do echo $i; done)

    for ((i=1; i <= $SAMPLES; i++)); do
        local output_file="${OUTDIR}/${OUTPUT_BASENAME}_${i}.tsv"

        if [[ -f ${output_file} && ! $FORCE == "yes" ]]; then
            echo "${output_file} exists, skipping!"
        else
            echo "Scaling tests: ${output_file}"
            for ii in $nodes; do
                echo "Running: $ii nodes"
                run_test $ii $output_file
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


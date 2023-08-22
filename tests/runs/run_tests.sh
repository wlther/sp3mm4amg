#!/bin/bash

settings_file="test_settings.conf"
results_dir="results"

# Read the first line to get the number of runs
num_runs=0
# Read each line from the file and execute programs
while read -r program; do
    if [[ $# -lt 5 ]]; then
        for threads in 1 2 4 8
        do
            export OMP_NUM_THREADS=$threads
            "./$program" $@
        done
    elif [[ $num_runs -gt 0 ]]; then
        echo $'\033[32;1m'"Implementation, R, AC, P, AC_NEXT, Preparation time, Threads, N-RAC, M-RAC, NNZ-RAC, Time-RAC, N-RACP, M-RACP, NNZ-RACP, Time-RACP"$'\033[0m'
        for threads in 1 2 4 8
        do
            export OMP_NUM_THREADS=$threads
            for (( i = 1; i <= num_runs; i++ )); do
            "./$program" $@ | tee -a "results.csv"
            done
        done
    else
        num_runs="$program"
    fi
done < "$settings_file"

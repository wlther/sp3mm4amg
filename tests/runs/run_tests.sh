#!/bin/bash

settings_file="test_settings.conf"
results_dir="results"
max_procs=$(nproc)

# Read the first line to get the number of runs
num_runs=0
# Read each line from the file and execute programs
while read -r program; do
    echo "$program"
    if [[ $# -lt 5 ]]; then
        # quick test with visual output in the terminal
        "./$program" $@
    elif [[ $num_runs -gt 0 ]]; then
        echo $'\033[32;1m'"Implementation,Collection,Size,Smoothing,Level,Preparation-time,Threads,N-1,M-1,NNZ-1,Time-1,N-2,M-2,NNZ-2,Time-2"$'\033[0m'
        threads=1
        while [[ $threads -lt $max_procs ]]
        do
            export OMP_NUM_THREADS=$threads
            export SP3MM_ITERATIONS=$num_runs
            "./$program" $@ | tee -a "$results_dir/results.csv"
            threads=$(( 2 * $threads ))
        done
    else
        num_runs="$program"
    fi
done < "$settings_file"

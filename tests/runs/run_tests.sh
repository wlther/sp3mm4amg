#!/bin/bash

settings_file="test_settings.conf"
results_dir="results"
suffix=$(date +'%d-%m-%Y-%H-%M-%S')

# Read the first line to get the number of runs
num_runs=0

# Read each line from the file and execute programs
while read -r program; do
    if [[ $num_runs -gt 0 ]]; then
        for (( i = 1; i <= num_runs; i++ )); do
            echo $'\033[32;1m'"$program"$'\033[0m'
            "./$program" $@ | tee "$results_dir/$program-$i"
        done
    else
        num_runs="$program"
    fi
done < "$settings_file"

#!/bin/bash
echo "Implementation,Collection,Size,Smoothing,Level,Preparation-time,Threads,N-1,M-1,NNZ-1,Time-1,N-2,M-2,NNZ-2,Time-2" > "results/results.csv"

while read line; do
	matrices=("$line")
	./sp3mm_benchmark ${matrices[@]} | tee -a "results/results.csv"
done < data/matrixGroups.list
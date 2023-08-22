#!/bin/bash
echo "Implementation,Collection,Size,Smoothing,Level,Preparation-time,Threads,N-RAC,M-RAC,NNZ-RAC,Time-RAC,N-RACP,M-RACP,NNZ-RACP,Time-RACP" > "results.csv"

while read line; do
	matrixes=("$line")
	./run_tests.sh ${matrixes[@]}
done < data/matrixGroups.list
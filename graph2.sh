#!/bin/bash

echo "for generating Eureka Graph"

echo "./experiment -N 262143 -K 1000"
 ./experiment -N 524288 -K 1000 -B 2 -E 0.5 -L 12 -l 5 -r 4 -t 1e-6 -e 1e-8  


echo "./experiment -N 262143 -K 1000"
 ./experiment -N 1048576 -K 1000 -B 2 -E 0.5 -L 12 -l 4 -r 4 -t 1e-6 -e 1e-8  


echo "./experiment -N 262143 -K 1000"
 ./experiment -N 2097152 -K 1000 -B 2 -E 0.2 -L 15 -l 3 -r 2 -t 1e-6 -e 1e-8  


echo "./experiment -N 262143 -K 1000"
 ./experiment -N 4194304 -K 1000 -B 2 -E 1  -L 12 -l 5 -r 4 -t 1e-6 -e 1e-8  


echo "./experiment -N 262143 -K 1000"
 ./experiment -N 8388608 -K 1000 -B 2 -E 1  -L 12 -l 5 -r 4 -t 1e-6 -e 1e-8

echo "./experiment -N 262143 -K 1000"
 ./experiment -N 1677216 -K 1000 -B 2 -E 1  -L 12 -l 5 -r 4 -t 1e-6 -e 1e-8
echo "./experiment -N 262143 -K 1000"
 ./experiment -N 33554432 -K 1000 -B 2 -E 1  -L 12 -l 5 -r 4 -t 1e-6 -e 1e-8
echo "./experiment -N 262143 -K 1000"
 ./experiment -N 67108864 -K 1000 -B 2 -E 1  -L 12 -l 5 -r 4 -t 1e-6 -e 1e-8
echo "./experiment -N 262143 -K 1000"
 ./experiment -N 134217728 -K 1000 -B 2 -E 1  -L 12 -l 5 -r 4 -t 1e-6 -e 1e

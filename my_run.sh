export OMP_NUM_THREADS=1
numactl --membind=0 --cpunodebind=0 ./experiment -N 33554432 -K 1000 -B 4 -E 4 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8 -M 256 -m 6

#for j in 50 1000 10000
#do
#    echo "K = $j"
#    ./experiment -N 32768 -K $j -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#    ./experiment -N 262144 -K $j -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#    ./experiment -N 1048576 -K $j -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#    ./experiment -N 4194304 -K $j -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#    ./experiment -N 8388608 -K $j -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#    ./experiment -N 16777216 -K $j -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#    ./experiment -N 33554432 -K $j -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#    #./experiment -N 268435456 -K $j -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#    #./experiment -N 1073741824 -K $j -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#done
# 

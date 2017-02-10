#echo "Sequential version..."

#./experiment -N 16777216 -K 100 -B 4 -E 2 -M 128 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#./experiment -N 131072 -K 50 -B 1 -E 1 -M 8 -m 2 -L 12 -l 4 -r 3 -t 1e-6 -e 1e-8
#./experiment -N 65536 -K 40 -B 4 -E 2 -M 128 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8

#for i in 32768 65536 131072 262144 524288 1048576 2097152 4194304 8388608 16777216 33554432
#do
#    for j in 50 1000 10000
#    do
#        time ./experiment -N $i -K $j -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#    done
#done

export BLOCK_SIZE=512
#for i in 1 4 8 12 16
#do
#    export OMP_NUM_THREADS=$i
#    echo "OMP = $i"
#    for j in 50 1000 10000
#    do
#        time ./experiment -N 32768 -K $j -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#        time ./experiment -N 262144 -K $j -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#        time ./experiment -N 1048576 -K $j -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#        time ./experiment -N 4194304 -K $j -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#        time ./experiment -N 8388608 -K $j -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#        time ./experiment -N 16777216 -K $j -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#        time ./experiment -N 33554432 -K $j -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#        time ./experiment -N 268435456 -K $j -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#        time ./experiment -N 1073741824 -K $j -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#    done
#done

for i in 1 2 4 8
do
    echo "OMP = $i"
    export OMP_NUM_THREADS=$i
    numactl --membind=0 --cpubind=0 ./experiment -N 33554432 -K 1000 -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
    #./experiment -N 33554432 -K 1000 -B 4 -E 4 -L 10 -l 6 -r 2 -M 256 -m 4 -t 1e-8 -e 1e-8
  done

for i in 12 16
do
    echo "OMP = $i"
    export OMP_NUM_THREADS=$i
    ./experiment -N 33554432 -K 1000 -M 256 -B 4 -E 4 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
done

#export OMP_NUM_THREADS=16
#./experiment -N 33554432 -K 1000 -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#numactl --membind=0 --cpubind=0 ./experiment -N 33554432 -K 1000 -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8
#./experiment -N 33554432 -K 1000 -M 256 -B 4 -E 2 -m 6 -L 10 -l 6 -r 2 -t 1e-8 -e 1e-8

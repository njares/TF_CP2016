#!/bin/bash

PROG_1=./grafo_completo/a.out

t_max=1
r_max=1

THREADS=12

# Header
echo "THREADS TIME RUTAS "

# Body
for ((t=$THREADS; t>=1; --t)); do
#for ((t=1; t<=$THREADS; ++t)); do
	# samples per #thread
	T=`OMP_NUM_THREADS=$t numactl --interleave=all $PROG_1 $t_max $r_max`
	D_x=`tail ./grafo_completo/data/output/Delta -n 1 | cut -d " " -f1`
	echo "$t $T $D_x"
done
echo " "

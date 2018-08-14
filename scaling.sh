#!/bin/bash

PROG_1=./grafo_completo/a.out

t_max=24
r_max=15

THREADS=12

# Header
echo -n "NAME"
for ((t=1; t<=$THREADS; ++t)); do
	echo -n " $t-threads"
done
echo " "

# Body
echo -n "mulatona "

for ((t=1; t<=$THREADS; ++t)); do
	# samples per #thread
	T=`OMP_NUM_THREADS=$t numactl --interleave=all $PROG_1 $t_max $r_max`
	echo -n "$T"
	echo -n " "
done
echo " "

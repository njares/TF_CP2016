#!/bin/bash

PROG_1=./grafo_completo/a.out
PROG_2=./min/a.out

t_max=3
r_max=30

D_x=`OMP_NUM_THREADS=12 $PROG_1 $t_max $r_max`
OMP_NUM_THREADS=12 numactl --interleave=all $PROG_2 $D_x

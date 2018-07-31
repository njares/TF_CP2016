#!/bin/bash

PROG_1=./grafo_completo/a.out
PROG_2=./min/a.out

t_max=3
r_max=30

T1=`OMP_NUM_THREADS=$1 $PROG_1 $t_max $r_max`
D_x=`wc -l grafo_completo/data/output/Delta | cut -d " " -f1`
T2=`OMP_NUM_THREADS=$1 numactl --interleave=all $PROG_2 $D_x`

echo "Tiempo de busqueda de rutas=$T1 segundos"
echo "Tiempo de minimización=$T2 segundos"

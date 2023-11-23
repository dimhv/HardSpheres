#!/bin/sh

gcc -std=c99  -Ofast -fopenmp 3d_hs_parallel_main.c -lm -o 3d_hs_parallel_main
# gcc -std=c99  -O3 -fopenmp 3d_hs_parallel_main.c -lm -o 3d_hs_parallel_main
# gcc -fopenmp 3d_hs_parallel_main.c -lm -o 3d_hs_parallel_main
# gcc 3d_hs_parallel_main.c -lm -o 3d_hs_parallel_main


